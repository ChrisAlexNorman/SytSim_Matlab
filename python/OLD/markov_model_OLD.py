import numpy as np
import numexpr as ne
import json
import jsbeautifier
from scipy.integrate import odeint, cumtrapz
from scipy.optimize import root_scalar
import multiprocessing as mp
import time
from numba import jit

def load_json_data(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    return data

class MarkovModel:


    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


    def __str__(self):
        return self.name + f": {self.n_states}-state Markov chain"


    def __len__(self):
        return self.n_states


    def _parse_model(self, model_dict):
        # Model Requirements
        for key in ["rates", "parameters", "initial_condition"]:
            if key not in model_dict:
                raise KeyError(f"{key} must be provided for model specification.")
        rates_shape = np.shape(model_dict["rates"])
        if len(rates_shape) != 2 or rates_shape[0] != rates_shape[1]:
            raise ValueError(f"rates is not square. Shape is {rates_shape}.")
        if len(model_dict["rates"]) != len(model_dict["initial_condition"]):
            raise ValueError(f"Rate matrix is {len(model_dict['rates'])}x{len(model_dict['rates'])} but initial condition has length {len(model_dict['initial_condition'])}.")
        
        # Derived/assumed properties
        model_dict.setdefault("n_states", len(model_dict["rates"]))
        model_dict.setdefault(
            "parameter_values",
            {name: info["value"] for name, info in model_dict["parameters"].items()},
        )
        model_dict.setdefault("simulations", [])
        model_dict.setdefault("name", "Unnamed")
        model_dict.setdefault("state_names", ["S" + str(n) for n in range(0, model_dict["n_states"])])
        model_dict.setdefault("multiprocessing", True)

        return model_dict


    def _parse_simulation(self, simulation):
        # Simulation requirements
        if "stimuli" not in simulation and "timestamp" not in simulation:
            raise ValueError("Either timestamped 'stimuli' or 'timestamp' must be specified.")
        simulation.setdefault("mode", "deterministic")
        recognised_modes = ["deterministic", "stochastic"]
        if simulation["mode"].lower() not in recognised_modes:
            raise ValueError(f"mode '{simulation['mode']}' not recognised. Must be one of {recognised_modes}.")
        
        # Derived/assumed properties
        simulation.setdefault("stimuli", {})
        simulation.setdefault("timestamp", [])
        simulation["aligned_stimuli_values"], simulation["timestamp"] = self._align_and_unpack_stimuli_values(simulation["stimuli"], simulation["timestamp"])
        simulation["sparse_transitions"] = self.sparsify_transitions(
                simulation["timestamp"], simulation["aligned_stimuli_values"]
            )
        simulation.setdefault("initial_condition", self.initial_condition)
        simulation.setdefault("state_names", self.state_names)
        simulation.setdefault("record", self.state_names)
        simulation.setdefault("runtime_cap", np.inf)
        simulation.setdefault("runtime", time.time())
        simulation.setdefault("n_simulations", 1)
        simulation.setdefault("n_processes", 1)
        
        return simulation


    def _align_and_unpack_stimuli_values(self, stimuli, timestamp):
        # If timestamp not specified, take timestamp from shortest stimulus.
        if len(timestamp) == 0:
            timestamp = [np.inf]
            for stimulus in stimuli.values():
                if stimulus['timestamp'][-1] < timestamp[-1]:
                    timestamp = stimulus['timestamp']

        # Extract stimulus values from stimulus data, interpolated at timestamp.
        stimuli_values = {
            name: np.interp(timestamp, stimulus['timestamp'], stimulus['value'])
            for name, stimulus in stimuli.items()
        }

        return stimuli_values, timestamp


    def _ode_system(self, states, t, timestamp, stimuli_values):
        '''Return rate of change of state probabilities, dp/dt, at time t given stimuli values.'''
        stimuli_at_t = {name: np.interp(t, timestamp, stimulus) for  name, stimulus in stimuli_values.items()}
        q_matrix_at_t = np.array(
            [
                [
                    ne.evaluate(str(element), {**self.parameter_values, **stimuli_at_t}).item()
                    for element in row
                ]
                for row in self.rates
            ]
        )
        q_matrix_at_t[np.diag_indices(self.n_states)] = -np.sum(q_matrix_at_t, axis=1)
        return np.dot(states, q_matrix_at_t)


    def import_model(self, filename):
        data = self._parse_model(load_json_data(filename))
        for key, value in data.items():
            setattr(self, key, value)
        return self


    def export_model(self, filename, data=None):
        if not data: data = {key: value for key, value in vars(self).items()}
        opts = jsbeautifier.default_options()
        opts.indent_size = 2
        data_beaut = jsbeautifier.beautify(json.dumps(data), opts)
        with open(filename, 'w') as file:
            file.write(data_beaut)
        return self


    def clear_simulations(self):
        simulations = self.simulations
        self.simulations = []
        return simulations


    def sparsify_transitions(self, timestamp, stimuli_values):
        '''
        Return dictionary of non-zero rates and their integrals evaluated at each timestamp for each source / destination pair.
        '''
        sparse_transitions = {}

        for source_idx, source_exit_rate_expressions in enumerate(self.rates):
            source_name = self.state_names[source_idx]

            for dest_idx, rate_expression in enumerate(source_exit_rate_expressions):
                dest_name = self.state_names[dest_idx]

                # Evaluate rate expression as a string equation given parameter and stimuli values
                rate_timeseries = ne.evaluate(str(rate_expression), {**self.parameter_values, **stimuli_values})

                # Only include non-zero transition rates
                if np.all(rate_timeseries == 0): continue

                # If the rate is a constant then transform it into a timeseries between the start and end times
                if rate_timeseries.ndim == 0:
                    rate_timestamp = np.array([0, timestamp[-1]])
                    rate_timeseries = np.array([rate_timeseries, rate_timeseries])
                else:
                    rate_timestamp = timestamp
                
                # Numerically integrate the transition rate
                rate_integral = cumtrapz(rate_timeseries, rate_timestamp, initial=0)

                # Append transition data
                transition = {
                    'destination' : dest_name,
                    'timestamp' : rate_timestamp,
                    'rate' : rate_timeseries,
                    'rate_integral' : rate_integral
                }
                if not source_name in sparse_transitions: sparse_transitions[source_name] = []
                sparse_transitions[source_name].append(transition)

        return sparse_transitions


    def jumps_to_probs(self, jump_time_dict, n_sims, timestamp):
        raw_probability = {
            name: {
                'timestamp': np.sort(jump_time),
                'probability': np.array(range(1,len(jump_time)+1))/n_sims
            }
            for name, jump_time in jump_time_dict.items()
        }
        probability = {
            name: np.interp(timestamp, np.sort(jump_time), np.array(range(1,len(jump_time)+1))/n_sims)
            for name, jump_time in jump_time_dict.items()
        }
        return probability, raw_probability


    def simulate(self, **kwargs):
        simulation = self._parse_simulation(kwargs)

        # Deterministic (ODE) solution
        if simulation['mode'].lower() == 'deterministic':
            # Numerically solve ODE
            states = odeint(
                self._ode_system,
                self.initial_condition,
                simulation["timestamp"],
                args=(
                    simulation["timestamp"],
                    simulation["aligned_stimuli_values"],
                ),
            )
            # Return results requested in 'record'   
            states = np.transpose(states)
            probability = {
                self.state_names[idx]: state
                for idx, state in enumerate(states)
                if self.state_names[idx] in simulation["record"]
            }
            simulation["probability"] = probability

        # Stochastic (Gillespie) solution
        elif simulation["mode"].lower() == "stochastic":
            # Initialise results container
            simulation["jump_time"] = {key: np.array([]) for key in simulation["record"]}

            with mp.Pool(simulation["n_processes"]) as pool:
                # Send jobs to worker pool
                async_results = [
                    pool.apply_async(stochastically_simulate, args=(simulation,))
                    for _ in range(simulation["n_simulations"])
                ]
                # Collect event times
                for output in async_results:
                    sim_result = output.get()
                    for key, values in sim_result.items():
                        simulation["jump_time"][key] = np.append(
                            simulation["jump_time"][key], (values)
                        )
            # Convert jump times to probabilities
            (
                simulation["probability"], simulation["raw_probability"]
            ) = self.jumps_to_probs(
                simulation["jump_time"], simulation["n_simulations"], simulation["timestamp"]
                )

        simulation['runtime'] = time.time() - simulation['runtime']
        self.simulations.append(simulation)
        return simulation

def wait_time_root_func(wait_time, shifted_rate_integrals, rand_log):
    '''
    Return sum of rate integrals up to wait_time plus given log(rand).
    Assumes integrals shifted to current_time = 0
    '''
    func_sum = rand_log
    for integral in shifted_rate_integrals:
        func_sum += np.interp(wait_time, integral['timestamp'], integral['value']) - integral['value'][0]

    return func_sum


def bisection(f, a, b, args,tol=1e-6, max_iter=100):
    """
    Bisection method to find a root of the function f in the interval [a, b].

    Parameters:
        f (callable): The function for which to find the root.
        a (float): The left endpoint of the interval.
        b (float): The right endpoint of the interval.
        tol (float): Tolerance for the root approximation (default: 1e-6).
        max_iter (int): Maximum number of iterations (default: 100).

    Returns:
        float: Approximation of the root.
    """
    for i in range(max_iter):
        c = (a + b) / 2
        if np.abs(f(c,args[0],args[1])) < tol:
            return c
        elif f(a,args[0],args[1]) * f(c,args[0],args[1]) < 0:
            b = c
        else:
            a = c
    raise RuntimeError("Bisection method did not converge.")

def stochastically_simulate(simulation):
    sparse_transitions = simulation["sparse_transitions"]
    current_time, *_, end_time = simulation["timestamp"]
    current_state = [simulation["state_names"][index] for index, value in enumerate(simulation["initial_condition"]) if value != 0]
    jump_time_dict = {key: [] for key in simulation['record']}
    while current_time < end_time  and time.time() - simulation['runtime'] < simulation['runtime_cap']:

        for state in current_state:
            if state in simulation['record']:
                jump_time_dict[state].append(current_time)

        if state not in sparse_transitions: break
        feasible_transitions = [transition|{'source': state} for state in current_state for transition in sparse_transitions[state]]

        # Find next jump time
        shifted_rate_integrals = [{
            'timestamp' : np.concatenate((
                [0], transition['timestamp'][transition['timestamp']>current_time]-current_time
                )),
            'value': np.concatenate((
                [np.interp(current_time, transition['timestamp'], transition['rate_integral'])],
                transition['rate_integral'][transition['timestamp']>current_time]
                ))
            }for transition in feasible_transitions]
        rand_log = np.log(np.random.random())
        try:
            # wait_time=bisection(wait_time_root_func,0,end_time-current_time,args=(shifted_rate_integrals, rand_log))
            wait_time = root_scalar(wait_time_root_func, x0=0, bracket=[0, end_time-current_time], args=(shifted_rate_integrals, rand_log), method='brentq').root # or 'toms748' 
        except:
            # No transition occured before end time
            break
        current_time += wait_time
        
        # find next jump state.
        current_rates = [np.interp(current_time, transition['timestamp'], transition['rate']) for transition in feasible_transitions]
        next_transition = np.searchsorted(np.cumsum(current_rates)/ np.sum(current_rates), np.random.random())
        for index, state in enumerate(current_state):
            if state == feasible_transitions[next_transition]['source']:
                current_state[index] = feasible_transitions[next_transition]['destination']
        
    return jump_time_dict
        