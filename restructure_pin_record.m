function Results = restructure_pin_record(Results)
%% Restructure raw record of pin releases from stochastic simulations
%
% Syntax: Results = restructure_pin_record(Results)
%
% Inputs:
%    Results::struct    Raw results of stochastic simulations
%
% Outputs:
%    Results::struct    Simulation results with restructured pin record
%
% m-file Requirements:
% -none-
%
% Author: Chris Norman, University of Warwick
% Source: https://github.com/ChrisAlexNorman/SytSim

%------------- BEGIN CODE --------------
fprintf('\nRestructuring record of SNAREpin releases:\n')
if ~isfield(Results,'PinRec')
    fprintf('!!! Pin data not recorded !!!\n')
else
    PinV = Results.PinRec;
    Results = rmfield(Results,'PinRec');
    nSNAREs = Results.metaData.nSNAREs;
    nVesicleSites = Results.metaData.nVesicleSites;
    
    % Initialise event times
    num_event_times = 0;
    for v = 1:length(PinV)
        num_event_times = num_event_times + size(PinV{v},1) - 1;
    end
    event_times = zeros(1,num_event_times);
    
    % Estimate variable sizes from event times
    freePins = zeros(nSNAREs+2,ceil(length(event_times)/3)+1);
    priPins = zeros(nSNAREs+2,ceil(length(event_times)/3)+1);
    triPins = zeros(nSNAREs+2,ceil(length(event_times)/3)+1);
    clear('event_times')
    Pins_on_release = zeros(2,length(Results.releaseTimes));
    freePins(2,1) = nVesicleSites; priPins(2,1) = nVesicleSites; triPins(2,1) = nVesicleSites;
    fpn = 2; ppn = 2; tpn = 2; rpn = 1;
    fprintf('Unpacking PinRec data... \nOf %d : ',length(PinV))
    for v = 1:length(PinV)
        % Display progress
        if ~mod(v,10000) || v == length(PinV)
            if v > 10000
                for j=1:length(num2str(v-1))
                    fprintf('\b'); % delete previous counter display
                end
            end
            fprintf('%d',v)
        end
        if ~(size(PinV{v},1) == 1)
            for ti = 2:size(PinV{v},1)
                if PinV{v}(ti,4) == nSNAREs+1
                    % Vesicle released
                    Pins_on_release(1,rpn) = PinV{v}(ti,1);
                    Pins_on_release(2,rpn) = PinV{v}(ti-1,4);
                    freePins(1,fpn) = PinV{v}(ti,1);
                    freePins(PinV{v}(ti-1,4)+2,fpn) = -1;
                    priPins(1,ppn) = PinV{v}(ti,1);
                    priPins(PinV{v}(ti-1,2)+2,ppn) = -1;
                    triPins(1,tpn) = PinV{v}(ti,1);
                    triPins(PinV{v}(ti-1,3)+2,tpn) = -1;
                    rpn=rpn+1; fpn=fpn+1; ppn=ppn+1; tpn=tpn+1;
                    
                else
                    if PinV{v}(ti,2) ~= PinV{v}(ti-1,2)
                        % primary Pin released
                        priPins(1,ppn) = PinV{v}(ti,1);
                        priPins(PinV{v}(ti,2)+2,ppn) = 1;
                        if PinV{v}(ti-1,4) <= nSNAREs
                            priPins(PinV{v}(ti-1,2)+2,ppn) = -1;
                        end
                        ppn = ppn + 1;
                    end
                    
                    if PinV{v}(ti,3) ~= PinV{v}(ti-1,3)
                        % tripartite Pin released
                        triPins(1,tpn) = PinV{v}(ti,1);
                        triPins(PinV{v}(ti,3)+2,tpn) = 1;
                        if PinV{v}(ti-1,4) <= nSNAREs
                            triPins(PinV{v}(ti-1,3)+2,tpn) = -1;
                        end
                        tpn = tpn + 1;
                    end
                    
                    if PinV{v}(ti,4) ~= PinV{v}(ti-1,4)
                        % check if both Pins released
                        freePins(1,fpn) = PinV{v}(ti,1);
                        freePins(PinV{v}(ti,4)+2,fpn) = 1;
                        if PinV{v}(ti-1,4) <= nSNAREs
                            freePins(PinV{v}(ti-1,4)+2,fpn) = -1;
                        end
                        fpn = fpn+1;
                    end
                end
            end
        end
    end
    clear('PinV')
    
    fprintf('  Done.\nRe-structuring freePins...  ')
    freePins(:,fpn:end)=[];
    try
        freePins = sortrows(freePins',1)';
    catch
        fprintf('\nRan out of memory. Trying slower method...  ')
        [freePins(1,:),sortfpn] = sort(freePins(1,:));
        for n = 1:nSNAREs
            freePins(n+1,:) = freePins(n+1,sortfpn);
        end
        clear('sortfpn');
    end
    
    fprintf('  Done.\nRe-structuring priPins...  ')
    priPins(:,ppn:end) = [];
    try
        priPins = sortrows(priPins',1)';
    catch
        fprintf('\nRan out of memory. Trying slower method...  ')
        [priPins(1,:),sortppn] = sort(priPins(1,:));
        for n = 1:nSNAREs
            priPins(n+1,:) = priPins(n+1,sortppn);
        end
        clear('sortppn');
    end
    
    fprintf('  Done.\nRe-structuring triPins...  ')
    triPins(:,tpn:end) = [];
    try
        triPins = sortrows(triPins',1)';
    catch
        fprintf('\nRan out of memory. Trying slower method...  ')
        [triPins(1,:),sorttpn] = sort(triPins(1,:));
        for n = 1:nSNAREs
            triPins(n+1,:) = triPins(n+1,sorttpn);
        end
        clear('sorttpn');
    end
    
    fprintf('  Done.\nRe-structuring Pins_on_release...  ')
    try
        Pins_on_release = sortrows(Pins_on_release',1)';
    catch
        fprintf('\nRan out of memory. Trying slower method...  ')
        [Pins_on_release(1,:),sortrpn] = sort(Pins_on_release(1,:));
        Pins_on_release(2,:) = Pins_on_release(2,sortrpn);
        clear('sortrpn');
    end
    Pins_on_release = Pins_on_release(2,:);
    
    fprintf('  Done.\nSumming rows...  ')
    for n = 1:nSNAREs+1
        freePins(n+1,:) = cumsum(freePins(n+1,:));
        priPins(n+1,:) = cumsum(priPins(n+1,:));
        triPins(n+1,:) = cumsum(triPins(n+1,:));
    end
    
    fprintf('  Done.\nThinning...')
    freePins = freePins(:,floor(linspace(1,size(freePins,2),min(size(freePins,2),1e6))));
    priPins = priPins(:,floor(linspace(1,size(priPins,2),min(size(priPins,2),1e6))));
    triPins = triPins(:,floor(linspace(1,size(triPins,2),min(size(triPins,2),1e6))));
    
    Results.priPins = priPins;
    Results.triPins = triPins;
    Results.freePins = freePins;
    Results.Pins_on_release = Pins_on_release;
    fprintf('  Done.\n')
end

metaData = Results.metaData;
Results = rmfield(Results,'metaData');
Results.metaData = metaData;
fprintf('All done!\n')
%------------- END OF CODE --------------

