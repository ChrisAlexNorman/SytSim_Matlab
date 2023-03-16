function [releaseTimes,repleniTimes,nVSim,PinRec] = simulate_vesicle(CaTimeSeries,recordPins,nSNAREs,P,kon,koff,kin,kouts,fmodel,R,rmodel,nRelPerWorker,timeCap,seed)
%% Repeatedly simulates specified model until target number of releases achieved
%
% Syntax:  [releaseTimes,repleniTimes,nVSim,PinRec] = simulate_vesicle(CaTimeSeries,recordPins,nSNAREs,P,kon,koff,kin,kouts,fmodel,R,rmodel,nRelPerWorker,timeCap,seed)
%
% Inputs:
%    CaTimeSeries::dbl  (N x 2) array: col1(time/ms) ; col2([Ca2+]/uM)
%    recordPins:bool    Toggle to record pin releases or not
%    nSNAREs::int       Number of SNARE complexes associated with vesicle
%    P::dbl             Prob. of any tripartite Syt being of type 1
%    kon::dbl           1/uMms; Ca-Syt binding rate
%    koff::dbl          1/ms;   Ca-Syt unbinding rate
%    kin::dbl           1/ms;   Syt-membrane association rate
%    kouts::dbl         1/ms;   Syt-membrane dissociation rates: (1)Syt1, (2)Syt7
%    fmodel::str        Fusion model name ("instant", "step" or "exponential")
%    R::dbl             Fusion rate associated parameter
%    rmodel::str        Vesicle replenishment model name ("fixed","delay" or "none")
%    nRelPerWorker::int Requested number of vesicle releases
%    timeCap::dbl       mins;    Real run time limit for vesicle simulation code
%    seed::int/str      rng seed
%
% Outputs:
%    releaseTimes::dbl - (1 x No.Releases) array of vesicle release times, ms
%    repleniTimes::dbl - (1 x No.Replens) array of vesicle replenishment times, ms
%    nVSim::int        - No. of vesicles simulated to achieve nRelPerWorker
%    PinRec::cell(dbl) - Record of number of released SNAREpins for each
%                        simulated vesicle as: [time, Pri, Tri, Both]
%
% Other m-files required:
% -none-
%
% See also: run_simulation
%
% Author: Chris Norman, University of Warwick
% Source: https://github.com/ChrisAlexNorman/SytSim

%------------- BEGIN CODE --------------
rng(seed)

%% Ca trace pre-calculations to reduce interpolation costs later
% Convert discrete [Ca2+] trace into piecewise polynomial object
CaTrace = spline(CaTimeSeries(:,1),CaTimeSeries(:,2)); %uM

% Numerically integrate discrete [Ca2+] time series
% and convert to piecewise polynomial object
CaIntDiscrete = cumtrapz(CaTimeSeries(:,1),CaTimeSeries(:,2));
CaInt = spline(CaTimeSeries(:,1),CaIntDiscrete);

% Extract useful values
TIME_END = CaTrace.breaks(end);
tSamples = linspace(0,TIME_END,1000);
CaIntSamples = ppval(CaInt,tSamples);

%%  Define vesicle for identicle initialisation based on inputs
% Define the type of the tripartite Syt at each SNARE
% 1 <-> Syt1; 2 <-> Syt7
if P > 1 % Fixed, non-probabilistic cases, always Syt7
    triSytType = 2*ones(1,nSNAREs);            
else
    triSytType = 2*ones(1,nSNAREs)-binornd(1,P,1,nSNAREs);
end

% Construct initial state and rate matrices
for s = nSNAREs:-1:1
    %% Primary Site
    if P == 2 % Primary site permanantly unoccupied (i.e. inserted)
        stateOccupancyCell{2*s-1}=[zeros(1,3),1];
        kout = 0;
    else
        stateOccupancyCell{2*s-1}=[1,zeros(1,3)];
        kout = kouts(1); % Always Syt1
    end
    rateConstantsCell{2*s-1} = [...
        0, 2*kon, 0, 0;...
        koff, 0, kon, 0;...
        0, 2*koff, 0, kin;...
        0, 0, kout, 0];
    CaDependentRatesCell{2*s-1} = [...
        0, 1, 0, 0;...
        0, 0, 1, 0;...
        0, 0, 0, 0;...
        0, 0, 0, 0];   
    %% Tripartite Site
    if P == 3 % Tripartite site permentantly unoccupied (i.e. inserted)
        stateOccupancyCell{2*s}=[zeros(1,3),1];
        kout = 0;
    else
        stateOccupancyCell{2*s}=[1,zeros(1,3)];
        kout = kouts(triSytType(s));
    end
    rateConstantsCell{2*s} = [...
        0, 2*kon, 0, 0;...
        koff, 0, kon, 0;...
        0, 2*koff, 0, kin;...
        0, 0, kout, 0];
    CaDependentRatesCell{2*s} = [...
        0, 1, 0, 0;...
        0, 0, 1, 0;...
        0, 0, 0, 0;...
        0, 0, 0, 0];
end

% Combine cells into arrays describing the vesicle for initialisation
statesINIT         = cell2mat(stateOccupancyCell);
rateConstantsFULL  = blkdiag(rateConstantsCell{:});
rateConstantsINIT  = sparse(rateConstantsFULL);
CaDependenciesINIT = sparse(blkdiag(CaDependentRatesCell{:}));

% Indicator arrays for membrane inserted states
insertionIndicesPri = 4:8:length(statesINIT);
insertionIndicesTri = 8:8:length(statesINIT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Begin main simulation loop                                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

releaseTimes=[]; repleniTimes=[]; nVSim=0; PinRec = cell(1);
tStart = tic;

% Exit main loop when sufficient releases achieved or max run time exceeded
while length(releaseTimes) < nRelPerWorker && toc(tStart) < 60*timeCap
    nVSim = nVSim + 1;
    nextTime = 0;
    if recordPins; PinRec{nVSim} = [0,0,0,0]; end
    
    % Simulate a vesicle 'site' up to calcium trace end time
    while nextTime < TIME_END
        if nextTime > 0 % Record vesicle replenishment
            repleniTimes = [repleniTimes,nextTime];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Vesicle initialisation
        states = statesINIT;
        rateConstants = rateConstantsINIT;
        CaDependencies = CaDependenciesINIT;
        
        nextIndex=[];
        priReleases = zeros(1,nSNAREs);
        triReleases = zeros(1,nSNAREs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulate a single vesicle until release or the end of time
        while true
            %% Update state vectors and rate matrices
            % Update switched states
            pastIndex = intersect(find(rateConstantsFULL(:,nextIndex)),find(states));
            states(nextIndex) = 1-states(nextIndex);
            states(pastIndex) = 1-states(pastIndex);
            
            % Create vector of current transition rates
            rates = states*rateConstants;
            ratesCaDep = states*CaDependencies;
            
            % Find released pins
            for s=1:nSNAREs
                priReleases(s) = any(find(states)==insertionIndicesPri(s));
                triReleases(s) = any(find(states)==insertionIndicesTri(s));
            end
            ReleasedPins  = priReleases.*triReleases;
            npriReleases = sum(priReleases);
            ntriReleases = sum(triReleases);
            nReleasedPins = sum(ReleasedPins);
            
            % Update pin record if the state of released pins has changed
            if recordPins
                if (npriReleases~=PinRec{nVSim}(end,2)) ||...
                   (ntriReleases~=PinRec{nVSim}(end,3)) ||...
                   (nReleasedPins~=PinRec{nVSim}(end,4))
                    PinRec{nVSim} = [PinRec{nVSim};nextTime,npriReleases,ntriReleases,nReleasedPins];
                end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Fusion model            
            switch fmodel
                case 'instant'
                    % Break loop immediately if required no. of pins released
                    if nReleasedPins >= R
                        break
                    else
                        fusionRate = 0;
                    end
                case 'step'
                    % Fusion has non-zero rate only for R and above released
                    if nReleasedPins >= R
                        fusionRate = 10; %1/ms
                    else
                        fusionRate = 0; %1/ms
                    end
                case 'exponential'
                    % A continuous function from the Arrhenius equation
                    fusionRate = exp(21.5)*1e-3 * exp(nReleasedPins*4.5-26); %1/ms
                otherwise
                    error(['Vesicle fusion model not recognised',newline,...
                'model must be "instant, "step", or "exponential"']);
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Direct Gillespie algorithm
            rates = [rates,fusionRate];
            ratesCaDep = [ratesCaDep,0];
            WD = rates*ratesCaDep';  % sum of Ca-dependent rate constants
            WI = rates*~ratesCaDep'; % sum of Ca-independent rate constants
            
            % Find next reaction TIME
            u1 = log(rand);
            Cat0 = ppval(CaInt,nextTime);
            itApprox = find(CaIntSamples > (WD*Cat0-WI*(tSamples-nextTime)-u1)/WD,1);
            try 
                nextTime = nextTime + fzero(@(tau) ...
                    WI*tau + WD*(ppval(CaInt,nextTime+tau)-Cat0) + u1,...
                    [tSamples(itApprox-1)-nextTime,tSamples(itApprox)-nextTime]);
            catch
                % No reaction occurs before TIME_END! Exit loop
                nextTime = TIME_END + 1;
                break
            end
            
            % Find next reaction INDEX
            Rcum = cumsum(rates*ppval(CaTrace,nextTime).*ratesCaDep + rates.*~ratesCaDep);
            nextIndex = find(Rcum >= rand*Rcum(end),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Exit Conditions
            if nextTime > TIME_END
                % Run time exceeded! Exit vesicle loop
                break
            elseif nextIndex == length(states)+1
                % Fusion occured! Exit vesicle loop
                break
            end
        end % of vesicle loop

        %% Store fusion time
        if nextTime <= TIME_END
            releaseTimes = [releaseTimes,nextTime];
            if recordPins
                PinRec{nVSim} = [PinRec{nVSim}; nextTime, npriReleases, ntriReleases, nSNAREs+1];
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Vesicle replenishment model
        switch rmodel
            case 'fixed'
                % Replenish after a fixed refractory time
                replenishmentTime = 2.5; %ms
            case 'delay'
                % fixed refractory time + random repriming time
                tRefractory = 2.5;  %ms refractory time
                krep = 0.02;        %1/ms repriming rate
                tRepriming  = exprnd(1/krep);
                replenishmentTime = tRefractory + tRepriming;
            case 'none'
                % Immediately exit loop
                break
            otherwise
                error(['Vesicle replenishment model not recognised',newline,...
                'model must be "fixed", "delay", or "none"']); 
        end
        nextTime = nextTime + replenishmentTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % of site loop
end % of main simulation loop
%------------- END OF CODE --------------
