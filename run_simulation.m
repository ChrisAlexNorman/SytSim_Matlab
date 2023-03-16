function Results = run_simulation(varargin)
%% Stochastically simulates synaptic vesicles for given scenario
%
% Syntax: Results = run_simulation(varargin)
%
% Inputs:
%    varargin::cell    Name-value pairs of parameter options
%
% Outputs:
%    Results::struct    See README.md
%
% m-file Requirements:
% simulate_vesicle.m
%
% See also: simulate_vesicle
%
% Author: Chris Norman, University of Warwick
% Source: https://github.com/ChrisAlexNorman/SytSim

%------------- BEGIN CODE --------------

%% Try to get or set up parallel pool
fprintf('Getting parallel pool...\n')

poolObj = gcp('nocreate');
if isempty(poolObj) % No pool already set up
    clear('poolObj');
    try % to set up pool with all available local workers
        myCluster = parcluster('local');
        nWorkers  = myCluster.NumWorkers;
        poolObj   = parpool('local',nWorkers);
    catch
        nWorkers  = 1;
        warning('Failed to start parallel pool')
    end
else
    nWorkers = poolObj.NumWorkers;
end

%% Declare default simulation parameters

% Parse inputs
if mod(nargin,2) ~= 0
    error('Name-value pair inputs only')
else
    inRequests = cell(2,nargin/2);
    i = 1;
    for iin = 1:2:nargin
        inRequests{1,i} = varargin{iin};
        inRequests{2,i} = varargin{iin+1};
        i = i+1;
    end
end

% Calcium trace as (N x 2) array: col1(time/ms) ; col2([Ca2+]/uM)
for i = 1:size(inRequests,2)
    if strcmp('CaTimeSeries',inRequests{1,i})
        CaTimeSeries = inRequests{2,i};
    end
end
if ~exist('CaTimeSeries','var')
    error('Ca2+ trace not provided in name-value format as: "CaTimeSeries", array')
end

fprintf('\nSimulation initialised with the following parameters:\n')

% rng seed
for i = 1:size(inRequests,2)
    if strcmp('seed',inRequests{1,i})
        seed = inRequests{2,i};
    end
end
if ~exist('Seed','var')
    seed = 'shuffle';
end
fprintf(['rng seed = ',num2str(seed),'\n'])

% Vesicle simulation run time cap (hours)
for i = 1:size(inRequests,2)
    if strcmp('timeCap',inRequests{1,i})
        timeCap = inRequests{2,i};
    end
end
if ~exist('timeCap','var')
    timeCap = 5; % mins
end
fprintf(['timeCap = ',num2str(timeCap),' mins\n'])

% Number of vesicle releases to collect
for i = 1:size(inRequests,2)
    if strcmp('nRelTot',inRequests{1,i})
        nRelTot = inRequests{2,i};
    end
end
if ~exist('nRelTot','var')
    nRelTot = 100000;
end
fprintf(['nRelTot = ',num2str(nRelTot),'\n'])

% Record pin releases?
for i = 1:size(inRequests,2)
    if strcmp('recordPins',inRequests{1,i})
        recordPins = inRequests{2,i};
    end
end
if ~exist('recordPins','var')
    recordPins = false;
end
fprintf(['recordPins = ',mat2str(recordPins),'\n'])

% Number of SNAREpins
for i = 1:size(inRequests,2)
    if strcmp('nSNAREs',inRequests{1,i})
        nSNAREs = inRequests{2,i};
    end
end
if ~exist('nSNAREs','var')
    nSNAREs = 6;
end
fprintf(['nSNAREs = ',num2str(nSNAREs),'\n'])

% Clamping model identifier
for i = 1:size(inRequests,2)
    if strcmp('P',inRequests{1,i})
        P = inRequests{2,i};
    end
end
if ~exist('P','var')
    P = 0; % Within [0,1] or equal to 2, 3, or 4
end
fprintf(['P = ',num2str(P),'\n'])

% Ca-binding 'on' rates (1/uMms)
for i = 1:size(inRequests,2)
    if strcmp('kon',inRequests{1,i})
        kon = inRequests{2,i};
    end
end
if ~exist('kon','var')
    kon = 1;
end
fprintf(['kon = ',num2str(kon),' 1/uMms\n'])

% Membrane insertion 'in' rates (1/ms)
for i = 1:size(inRequests,2)
    if strcmp('kin',inRequests{1,i})
        kin = inRequests{2,i};
    end
end
if ~exist('kin','var')
    kin = 100;
end
fprintf(['kin = ',num2str(kin),' 1/ms\n'])

% Vesicle fusion model name
for i = 1:size(inRequests,2)
    if strcmp('fmodel',inRequests{1,i})
        fmodel = inRequests{2,i};
    end
end
if ~exist('fmodel','var')
    fmodel = 'exponential'; % from 'instant','step','exponential'
end
fprintf(['fmodel = ',fmodel,'\n'])

% Vesicle fusion related parameter
for i = 1:size(inRequests,2)
    if strcmp('R',inRequests{1,i})
        R = inRequests{2,i};
    end
end
if ~exist('R','var')
    R = 3; % 0 to nSNAREs
end
fprintf(['R = ',num2str(R),'\n'])

% Vesicle replenishment model name
for i = 1:size(inRequests,2)
    if strcmp('rmodel',inRequests{1,i})
        rmodel = inRequests{2,i};
    end
end
if ~exist('rmodel','var')
    rmodel = 'none'; % from 'fixed','delay','none'
end
fprintf(['rmodel = ',rmodel,'\n'])

%% find koff and kouts for Syt1 and Syt7

koff = kon * 150;
fprintf(['koff = ',num2str(koff),' 1/ms\n'])

kouts(1) = 0.5*(0.5-kin-2*koff)/(0.5-2*koff);
fprintf(['Syt1 kout = ',num2str(kouts(1)),' 1/ms\n'])

kouts(2) = 0.015*(0.015-kin-2*koff)/(0.015-2*koff);
fprintf(['Syt7 kout = ',num2str(kouts(2)),' 1/ms\n'])

%% Evaluate Results

% Number of releases to be collected by each worker
nRelPerWorker = ceil(nRelTot / nWorkers);

% Initialise results objects
releaseTimes = cell(1,nWorkers);
repleniTimes = cell(1,nWorkers);
nVsim = zeros(1,nWorkers);
PinRecRaw = cell(1,nWorkers);
evaluationFutures(nWorkers) = parallel.FevalFuture();

fprintf('\nSending jobs to workers...\n')

% Send tasks to parallel workers
tStart = tic;
for v = 1:nWorkers
    evaluationFutures(v) = parfeval(poolObj,@simulate_vesicle,4,CaTimeSeries,recordPins,nSNAREs,P,kon,koff,kin,kouts,fmodel,R,rmodel,nRelPerWorker,timeCap,seed);
    pause(0.1) % for rng based on system clock
end

fprintf('Done. Waiting for results...\n')

% Collect results from workers as they arrive
for vi = 1:nWorkers
    [~,releaseTimes{vi},repleniTimes{vi},nVsim(vi),PinRecRaw{vi}] = fetchNext(evaluationFutures);
end
tElapsed = toc(tStart) / 60;

fprintf(['Done! Elapsed time = ',num2str(tElapsed),' mins\n'])

%% Construct output
fprintf('\nAggregating results... \n')

% Aggregate results from workers
releaseTimes = cell2mat(releaseTimes);
repleniTimes = cell2mat(repleniTimes);
nVesicleSites = sum(nVsim);
PinRec = cell(0);
for w = 1:nWorkers
    PinRec = [PinRec,PinRecRaw{w}];
end

% Collect simulation info into a metaData structure
metaData.CaTimeSeries = CaTimeSeries;
metaData.nVesicleSites = nVesicleSites;
metaData.nSNAREs = nSNAREs;
metaData.P = P;
metaData.kon = kon;
metaData.koff = koff;
metaData.kin = kin;
metaData.kouts = kouts;
metaData.fmodel = fmodel;
metaData.R = R;
metaData.rmodel = rmodel;
metaData.seed = seed;
metaData.nWorkers = nWorkers;
metaData.runTime = tElapsed;
metaData.timeCap = timeCap;

% Build output Results structure
Results.releaseTimes = releaseTimes;
Results.repleniTimes = repleniTimes;
if recordPins; Results.PinRec = PinRec; end
Results.metaData = metaData;

% Close parallel pool
if exist('poolObj','var')
    delete(poolObj);
end
fprintf('Done.\n')

fprintf('\nSimulations complete!\n')

end
%------------- END OF CODE --------------
