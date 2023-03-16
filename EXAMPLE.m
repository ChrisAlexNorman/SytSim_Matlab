%% Example script to extract responses of several models to an AP pair
%
% m-file Requirements:
% run_simulation.m
% simulate_vesicle.m
%
% See also: run_simulation, simulate_vesicle
%
% Author: Chris Norman, University of Warwick
% Source: https://github.com/ChrisAlexNorman/SytSim

% Directory and file prefix to save results under
SavePath = './Results/2AP_40nm_20ms';

% Load Ca2+ trace
CaTimeSeries = csvread('./CalciumTraces/2AP_40nm_20ms.csv');

% Parameters to vary
Ps = [0, 1, 3]; % Three limiting cases of clamping architecture
rmodels = {'none','delay'}; % With and without vesicle replenishment

%% Run and save simulations

for P = Ps
    for rmodel_idx = 1:length(rmodels)
        rmodel = rmodels{rmodel_idx};
        fprintf(['---------- STARTING: P = ',num2str(P),', rmodel = ', rmodel,' ----------\n\n'])

        % Collect raw Results
        Results = run_simulation('CaTimeSeries', CaTimeSeries, 'P', P, 'rmodel', rmodel, 'recordPins', true);

        % Restructure record of SNAREpin releases in Results
        Results = restructure_pin_record(Results);

        % Save results
        save([SavePath,'_P-',num2str(P),'_rmodel-',rmodel,'.mat'],'-struct','Results');
        
        fprintf(['\n---------- FINISHED: P = ',num2str(P),', rmodel = ', rmodel,' ----------\n\n'])
    end
end
