% Nguyen et al., "Coinfecting Phages Impede Each Other's Entry Into The Cell", Current Biology 2024
%% 1. PREPARATION OF PARAMETERS AND DATA STRUCTURE
clc
close all
clear all

rng('shuffle');

%Simulation setup
nVector = 0:20; %The values of outside MOI to be simulated
numCellsPerMOI = 10000;
tAds = 1 * 60; %seconds, when phages adsorb to the cell with respect to phage perfusion; not really used in the analysis later

%Function for drawing the waiting time later, as in the Gillespie algorithm
func_inverseCDF = @(a, r) 1./a .* log(1./(1 - r)) .* (r >= 0) .* (r <= 1);

%Function to parametrize the kinetic parameters (eta, k, tau)
func_expDecay = @(a, b, c, x) a .* exp(b.*(x-1)) + c;

%Kinetic parameters (MOI-dependent or not), obtained from microfluidic infection data
paramFolder = [pwd '/Dependent files/'];
paramFileName = 'paramInferred-14-Mar-2023.mat';
load([paramFolder paramFileName]);

flagMOIdependence = 1; %Set by user; by default, MOI-dependent parameterization
if flagMOIdependence == 0 %MOI-independent
    func_eta = @(n) 0.50;
    func_k = @(n) 0.01;
    func_tau = @(n) 30;
elseif flagMOIdependence == 1 %Parametrized smoothly, MOI is the only remaining input
    func_eta = @(n) func_expDecay(paramInferred.eta.param(1), paramInferred.eta.param(2), paramInferred.eta.param(3), n);
    func_k = @(n) func_expDecay(paramInferred.k.param(1), paramInferred.k.param(2), paramInferred.k.param(3), n);
    func_tau = @(n) paramInferred.tau.param;
end

%Initialization of data structure for all cells
simTimeVector = 0:1:(60*30); %seconds
simData = struct();

%Instruction for plotting
modelColor = [234 70 71]/255;
approxColor = [38 79 209]/255;

disp('Preparation for the simulation is complete.');

%% 2. SIMULATION AND SAVING RESULTS (No plotting)
clc
close all

timerStart = tic;
ii = 1; %Master cell index

%Simulation for each outside MOI
fprintf('Now simulating for outside MOI =');
for i_n = 1:numel(nVector)
    %Specification and display
    n = nVector(i_n);
    fprintf(' %d', n);

    %Kinetic parameters for this outside MOI, shared by all cells
    eta = func_eta(n);
    k = func_k(n);
    tau = func_tau(n);

    %Simulation of the time-series for each cell
    for i_cell = 1:numCellsPerMOI
        %Number of injectable phages among those that adsorb
        m = binornd(n, eta);

        %Initialization
        outside = m;
        inside = 0;
        currentTime = 0;
        vectorTimeBetweenInitiations = zeros(1, m); %One element for the entry time of each phage

        %Gillespie algorithm, for the waiting time between phage entries
        while inside < m %outside > 0
            %Propensities of phage entries among the adsorbed phages that have not entered the cell yet
            propensity = (m - inside)*k; %outside * k

            %Generating the waiting time to entry initiation, from a random number ~Unif[0,1]
            waitingTime = func_inverseCDF(propensity, rand(1, 1));

            %Saving time information
            vectorTimeBetweenInitiations(inside+1) = waitingTime;

            %Updating the time and counts of phages
            currentTime = currentTime + waitingTime;
            outside = outside - 1;
            inside = inside + 1;
        end

        %Saving information for this cell, in the same organization as cellData from experiments
        simData(ii).exp = 'Simulation';
        simData(ii).timeResolution = simTimeVector(2)-simTimeVector(1); %Imaging resolution, here, as if imaged very fast
        simData(ii).outsideMaxCount = n;
        simData(ii).insideMaxCount = m; %All injectable phages would enter the cell eventually

        %With respect to the time of phage adsorption
        %Waiting time until cell entry initiations for this cell
        simData(ii).waitingTimeBtwn = vectorTimeBetweenInitiations; %Not exist in cellData
        simData(ii).waitingTimeSinceAds = cumsum(vectorTimeBetweenInitiations); %Not exist in cellData

        %Entry time since phage adsorption
        simData(ii).entryTimesAfterAds = simData(ii).waitingTimeSinceAds + tau;

        %Intracellular phage number, with respect to time since phage adsorption
        simData(ii).biolTimeVector = simTimeVector;
        simData(ii).insideVectorBiolTime = zeros(1, numel(simData(ii).biolTimeVector));
        for i_entry = 1:numel(simData(ii).entryTimesAfterAds)
            flagTime = (simData(ii).biolTimeVector >= simData(ii).entryTimesAfterAds(i_entry));
            simData(ii).insideVectorBiolTime(flagTime) = i_entry;
        end

        %With respect to the time of microfluidic perfusion
        %Entry time since microfluidic perfusion
        simData(ii).injPhageTime = simData(ii).entryTimesAfterAds + tAds;

        %Intracellular phage number, with respect to time since microfluidic perfusion
        simData(ii).flowTimeVector = simTimeVector;
        simData(ii).insideVectorFlowTime = zeros(1, numel(simData(ii).flowTimeVector));
        for i_entry = 1:numel(simData(ii).injPhageTime)
            flagTime = (simData(ii).flowTimeVector >= simData(ii).injPhageTime(i_entry));
            simData(ii).insideVectorFlowTime(flagTime) = i_entry;
        end

        %Finally, setting up adsorptions to be synchronous
        simData(ii).adsPhageTime = repelem(tAds, 1, simData(ii).outsideMaxCount);
        simData(ii).adsTimeRange = range(simData(ii).adsPhageTime); %0 by setup
        simData(ii).outsideVectorFlowTime = zeros(1, numel(simData(ii).flowTimeVector));
        flagTime = (simData(ii).flowTimeVector >= tAds);
        simData(ii).outsideVectorFlowTime(flagTime) = tAds;

        %Move forward in the master cell index
        ii = ii+1;
    end
end

timerEnd = toc(timerStart);

disp([newline 'Simulation and calculation are complete, taking ' num2str(timerEnd, '%.2f') ' sec.']);

%% 3.1 TIME-DEPENDENT INSIDE MOI FOR A GIVEN OUTSIDE MOI
clc
close all
clear graphSim graphAna

%For this observable, the analytical expression is simple enough that there is no need for any phenomenological approximation.

%Further preparations
nRow = 3;
nCol = 4;
outsideToShow = 1:12;

%For simulation data, show as discrete markers
timeSimToShow = (0:20)*60; %Seconds
flagTimeSimToShow = ismember(simTimeVector, timeSimToShow);

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_outside = 1:numel(outsideToShow)
    subplot(nRow, nCol, i_outside);
    hold on;

    currentOutside = outsideToShow(i_outside);
    
    %SIMULATION
    flagOutside = [simData(:).outsideMaxCount] == currentOutside;
    currentData = simData(flagOutside);

    %Preparing the stacked array, and calculating statistics
    stackedInsideMOIOverTime = zeros(numel(currentData), numel(simTimeVector));
    for i_cell = 1:numel(currentData)
        stackedInsideMOIOverTime(i_cell, :) = [currentData(i_cell).insideVectorBiolTime];
    end
    statsInsideMOIOverTime(1, :) = mean(stackedInsideMOIOverTime, 1);

    %ANALYTICAL
    anaInsideTimeVector = zeros(1, numel(simTimeVector));
    inputs = [repelem(currentOutside, numel(simTimeVector), 1), simTimeVector'];
    params = [func_eta(currentOutside) func_k(currentOutside) func_tau(currentOutside)];
    anaInsideTimeVector = funcInsideOutside_v2(params, inputs);

    %VISUALIZATION
    graphAna = plot(simTimeVector/60, anaInsideTimeVector,...
        'Marker', 'None', 'LineStyle', '-', 'LineWidth', 1.5, 'Color', modelColor,...
        'DisplayName', ['Analytical predictions']);

    graphStats = plot(simTimeVector(flagTimeSimToShow)/60, statsInsideMOIOverTime(1, flagTimeSimToShow),...
        'LineStyle', 'None', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'DisplayName', ['Simulation mean']);

    %Adjusting plot properties
    if currentOutside == 1
        legend([graphStats graphAna], 'FontSize', 8, 'Location', 'SouthEast');
    end
    pbaspect([1.62 1 1])
    xlim([0 10.5]);
    ylim([0 currentOutside]);
    xlabel('Time (min)');
    ylabel('Inside MOI');
    title(['Outside MOI = ' num2str(currentOutside)], 'FontSize', 12);
    grid on; box on;
end

%% 3.2 SCALING BETWEEN INSIDE AND OUTSIDE MOI AT A GIVEN TIME
clc
close all
clear graphSim graphAna

%Approximation fit
xspace = 0:0.1:30; %in time, minutes
approxFit.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf], 'Upper', [Inf Inf], 'Startpoint', [1 1]);
approxFit.function = @(a, K, n) a .* n ./ (K + n);
approxFit.model = fittype(approxFit.function, 'Options', approxFit.options,...
    'Independent', 'n', 'Coefficients', {'a', 'K'});

%Further preparations
nRow = 3;
nCol = 4;
timeInterest = 0:1:10; %minutes
outsideToShow = 0:12;

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_time = 1:numel(timeInterest)
    subplot(nRow, nCol, i_time);
    hold on;

    currentTime = timeInterest(i_time);
    flagTime = (simTimeVector == currentTime*60);

    %SIMULATION
    inoutPairs = zeros(numel(simData), 2);
    inoutPairs(:, 1) = [simData(:).outsideMaxCount];
    for i_cell = 1:numel(simData)
        inoutPairs(i_cell, 2) = [simData(i_cell).insideVectorBiolTime(flagTime)];
    end
    insideAverage = zeros(1, numel(outsideToShow));
    for i_outside = 1:numel(outsideToShow)
        currentOutside = outsideToShow(i_outside);
        flagOutside = inoutPairs(:, 1) == currentOutside;

        insideAverage(i_outside) = mean(inoutPairs(flagOutside, 2));
    end

    %ANALYTICAL
    anaInsideVector = zeros(1, numel(outsideToShow));
    for i_outside = 1:numel(outsideToShow)
        currentOutside = outsideToShow(i_outside);
        inputs = [currentOutside currentTime*60]; %Model in seconds
        params = [func_eta(currentOutside) func_k(currentOutside) func_tau(currentOutside)];
        anaInsideVector(i_outside) = funcInsideOutside_v2(params, inputs);
    end

    %APPROXIMATION
    [xToFit, yToFit] = prepareCurveData(outsideToShow, insideAverage);
    [param, ~] = fit(xToFit, yToFit, approxFit.model); %Fitting in minutes
    yspace = approxFit.function(param.a, param.K, xspace);

    %VISUALIZATION
    graphAna = plot(outsideToShow, anaInsideVector,...
        'Marker', 'None', 'LineStyle', '-', 'LineWidth', 1.5, 'Color', modelColor,...
        'DisplayName', ['Analytical predictions']);
 
    graphApprox = plot(xspace, yspace,...
        'Marker', 'None', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', approxColor,...
        'DisplayName', ['Approximation']);

    graphStats = plot(outsideToShow, insideAverage,...
        'LineStyle', 'None', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'DisplayName', ['Simulation mean']);

    %Adjusting plot properties
    if i_time == 1
        legend([graphStats graphAna graphApprox], 'FontSize', 8, 'Location', 'NorthWest');
    end
    pbaspect([1 1 1])
    xlim([0 12]);
    ylim([0 5]);
    xlabel('Outside MOI');
    ylabel('Inside MOI');
    title(['At \it\rm\bf = ' num2str(currentTime) ' min'], 'FontSize', 12);
    grid on; box on;
end

%% 3.3 ENTRY EFFICIENCY AS A FUNCTION OF OUTSIDE MOI
clc
close all
clear graphSim graphAna

%Approximation fit
xspace = 0:0.1:30; %in time, minutes
approxFit.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [-Inf -Inf], 'Upper', [Inf Inf], 'Startpoint', [1 1]);
approxFit.function = @(a, K, n) a .* n ./ (K + n);
approxFit.model = fittype(approxFit.function, 'Options', approxFit.options,...
    'Independent', 'n', 'Coefficients', {'a', 'K'});

%Further preparations
nRow = 2;
nCol = 5;
timeInterest = 1:1:10; %minutes
outsideToShow = 0:12;

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_time = 1:numel(timeInterest)
    subplot(nRow, nCol, i_time);
    hold on;

    currentTime = timeInterest(i_time);
    flagTime = (simTimeVector == currentTime*60);

    %SIMULATION
    inoutPairs = zeros(numel(simData), 2);
    inoutPairs(:, 1) = [simData(:).outsideMaxCount];
    for i_cell = 1:numel(simData)
        inoutPairs(i_cell, 2) = [simData(i_cell).insideVectorBiolTime(flagTime)];
    end
    insideAverage = zeros(1, numel(outsideToShow));
    for i_outside = 1:numel(outsideToShow)
        currentOutside = outsideToShow(i_outside);
        flagOutside = inoutPairs(:, 1) == currentOutside;

        insideAverage(i_outside) = mean(inoutPairs(flagOutside, 2));
    end

    %ANALYTICAL
    anaInsideVector = zeros(1, numel(outsideToShow));
    for i_outside = 1:numel(outsideToShow)
        currentOutside = outsideToShow(i_outside);
        inputs = [currentOutside currentTime*60]; %Model in seconds
        params = [func_eta(currentOutside) func_k(currentOutside) func_tau(currentOutside)];
        anaInsideVector(i_outside) = funcInsideOutside_v2(params, inputs);
    end

    %APPROXIMATION
    [xToFit, yToFit] = prepareCurveData(outsideToShow, insideAverage);
    [param, ~] = fit(xToFit, yToFit, approxFit.model); %Fitting in minutes
    yspace = approxFit.function(param.a, param.K, xspace);

    %VISUALIZATION
    graphAna = plot(outsideToShow, anaInsideVector./outsideToShow,...
        'Marker', 'None', 'LineStyle', '-', 'LineWidth', 1.5, 'Color', modelColor,...
        'DisplayName', ['Analytical predictions']);
 
    graphApprox = plot(xspace, yspace./xspace,...
        'Marker', 'None', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', approxColor,...
        'DisplayName', ['Approximation']);

    graphStats = plot(outsideToShow, insideAverage./outsideToShow,...
        'LineStyle', 'None', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None',...
        'DisplayName', ['Simulation mean']);

    %Adjusting plot properties
    if i_time == 1
        legend([graphStats graphAna graphApprox], 'FontSize', 8, 'Location', 'NorthWest');
    end
    pbaspect([1 1 1])
    xlim([0 12]);
    xticks(0:3:12);
    ylim([0 0.8]);
    yticks(0:0.2:0.8);
    xlabel('Outside MOI');
    ylabel('Entry efficiency');
    title(['At \it\rm\bf = ' num2str(currentTime) ' min'], 'FontSize', 12);
    grid on; box on;
end

%% 3.4 DISTRIBUTION OF INSIDE MOI GIVEN AN OUTSIDE MOI AT A GIVEN TIME
clc
close all
clear graphSim graphAna

%Approximation fit
approxFit.options = fitoptions('Method', 'NonlinearLeastSquares',...
    'Lower', [0], 'Upper', [Inf], 'Startpoint', [1]);
approxFit.function = @(mu, n, lambda) sum(poisspdf(0:1:n, mu))^(-1) * poisspdf(lambda, mu);
approxFit.model = fittype(approxFit.function, 'Options', approxFit.options,...
    'Independent', 'lambda', 'Coefficient', 'mu', 'Problem', 'n');

%Further preparations
flagShowGraph = {'Off' 'On'};

timeToShow = [1 2 5 10];
outsideToShow = [1 2 5 10 20];

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

for i_outside = 1:numel(outsideToShow)
    %MOI-specific preparations
    currentOutside = outsideToShow(i_outside);
    flagOutside = ([simData(:).outsideMaxCount] == currentOutside);
    insideValues = 0:1:currentOutside;
    binMOI = (-0.5):1:(currentOutside+0.5);
    xspace = 0:1:currentOutside;

    for i_time = 1:numel(timeToShow)
        subplot(numel(timeToShow), numel(outsideToShow), (i_time-1)*numel(outsideToShow) + i_outside);
        hold on;

        currentTime = timeToShow(i_time);

        %SIMULATION
        flagTime = (simTimeVector == currentTime*60);
        insideVector = [];
        for i_cell = find(flagOutside)
            insideVector = [insideVector simData(i_cell).insideVectorBiolTime(flagTime)];
        end
        [freq, edges] = histcounts(insideVector, binMOI, 'Normalization', 'PDF');
        
        %ANALYTICAL
        anaInsideProbVector = zeros(1, numel(insideValues));
        params = [func_eta(currentOutside) func_k(currentOutside) func_tau(currentOutside)];
        for i_inside = 1:numel(insideValues)
            currentInside = insideValues(i_inside);
            inputs = [currentOutside currentInside currentTime*60];
            anaInsideProbVector(i_inside) = funcProbInsideOutside_v2(params, inputs);
        end

        %APPROXIMATION
        [xToFit, yToFit] = prepareCurveData(insideValues, freq);
        [param, ~] = fit(xToFit, yToFit, approxFit.model, 'problem', currentOutside);
        yspace = approxFit.function(param.mu, currentOutside, xspace);

        %VISUALIZATION
        graphSim = bar(insideValues, freq, 0.8,...
            'FaceColor', 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'None',...
            'DisplayName', 'Simulation');

        graphAna = plot(insideValues, anaInsideProbVector,...
            'Marker', 'None', 'MarkerFaceColor',  modelColor, 'MarkerEdgeColor', 'None', 'MarkerSize', 4,...
            'LineStyle', '-', 'Color', modelColor, 'LineWidth', 1.5,...
            'DisplayName', 'Model');

        graphApprox = plot(xspace, yspace,...
            'Marker', 'None', 'LineStyle', '--', 'LineWidth', 1.5, 'Color', approxColor,...
            'DisplayName', ['Approximation']);

        %Adjusting plot properties
        if i_time == 1 & i_outside == 1
            legend([graphSim graphAna graphApprox], 'Location', 'NorthEast');
        end
        xlim([min(binMOI) max(binMOI)]);
        ylim(get(gca, 'YLim')*1.1);
        xlabel('Inside MOI');
        ylabel('Probability');
        title(['Outside MOI = ' num2str(currentOutside) ' @ ' num2str(currentTime) ' min'], 'FontSize', 10);
        grid on; box on;
    end
end

%% 3.5 PDF OF ENTRY TIME OF A GIVEN ENTRY ORDER FOR A GIVEN OUTSIDE MOI
clc
close all
clear graphSim graphAna

%Further preparations
outsideToShow = [1 2 5 10 20];
entriesToShow = [1:5];

figure('Units', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
hold on;

i_subplot = 1;
for i_entry = 1:numel(entriesToShow)
    currentEntry = entriesToShow(i_entry);

    for i_outside = 1:numel(outsideToShow)
        currentOutside = outsideToShow(i_outside);

        %Invalid if the inquired phage entry is more than the number of adsorbed phages
        if currentOutside < currentEntry
            i_subplot = i_subplot+1;
            continue
        end

        subplot(numel(entriesToShow), numel(outsideToShow), i_subplot);
        hold on;

        %SIMULATION
        flagOutside = [simData(:).outsideMaxCount] == currentOutside;
        flagInside = [simData(:).insideMaxCount] >= currentEntry; %Permiting this entry
        currentData = simData(flagOutside & flagInside);

        %Concatenate the entry time of this given entry order in these cells
        concatEntryTime = zeros(1, numel(currentData));
        for i_cell = 1:numel(currentData)
            concatEntryTime(i_cell) = currentData(i_cell).entryTimesAfterAds(currentEntry);
        end

        %ANALYTICAL
        anaPDFEntryTime = zeros(1, numel(simTimeVector));
        inputs = [repelem(currentOutside, numel(simTimeVector), 1),...
            repelem(currentEntry, numel(simTimeVector), 1),...
            simTimeVector'];
        params = [func_eta(currentOutside) func_k(currentOutside) func_tau(currentOutside)];
        anaPDFEntryTime = funcPDFEntryTime_v2(params, inputs);
        renormalization = trapz(simTimeVector/60, anaPDFEntryTime'); %Calculated in seconds, but displayed in minutes

        %Visualization
        graphSim = histogram(concatEntryTime/60, (0:10:max(concatEntryTime))/60, 'Normalization', 'PDF',...
            'FaceColor', 'k', 'EdgeColor', 'None', 'FaceAlpha', 0.25,...
            'DisplayName', 'Simulation');

        graphAna = plot(simTimeVector/60, anaPDFEntryTime'/renormalization,...
            'Marker', 'None', 'LineStyle', '-', 'LineWidth', 1, 'Color', 'r',...
            'DisplayName', 'Model');

        %Adjusting plot properties
        if i_subplot == 1
            legend([graphSim graphAna], 'Location', 'NorthEast', 'FontSize', 7);
        end
        pbaspect([1.62 1 1])
        xlabel('Entry time (min)');
        ylabel('Prob. density');
        title(['Entry #' num2str(currentEntry) ' for Outside MOI = ' num2str(currentOutside)]);
        xlim([0 12]);
        xticks(0:4:12);
        ylim(get(gca, 'YLim')*1.1);
        grid on; box on;

        i_subplot = i_subplot+1;
    end
end