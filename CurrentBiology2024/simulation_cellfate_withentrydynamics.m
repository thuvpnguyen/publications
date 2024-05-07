% Nguyen et al., "Coinfecting Phages Impede Each Other's Entry Into The Cell", Current Biology 2024
%% 1.1 SIMULATION OF PHAGE ENTRY DYNAMICS
clc
close all
clear all

rng('shuffle');

%Simulation setup
nVector = 0:5; %The values of outside MOI to be simulated
numCellsPerMOI = 1000;
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
biolTimeVector = 0:1:(60*30); %seconds
simData = struct();

disp('Preparation for the entry dynamics simulation is complete.');

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

        %Instantiation
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
        simData(ii).timeResolution = biolTimeVector(2)-biolTimeVector(1); %Imaging resolution, here, as if imaged very fast
        simData(ii).outsideMaxCount = n;
        simData(ii).insideMaxCount = m; %All injectable phages would enter the cell eventually

        %With respect to the time of phage adsorption
        %Waiting time until cell entry initiations for this cell
        simData(ii).waitingTimeBtwn = vectorTimeBetweenInitiations; %Not exist in cellData
        simData(ii).waitingTimeSinceAds = cumsum(vectorTimeBetweenInitiations); %Not exist in cellData

        %Entry time since phage adsorption
        simData(ii).entryTimesAfterAds = simData(ii).waitingTimeSinceAds + tau;

        %Intracellular phage number, with respect to time since phage adsorption
        simData(ii).biolTimeVector = biolTimeVector;
        simData(ii).insideVectorBiolTime = zeros(1, numel(simData(ii).biolTimeVector));
        for i_entry = 1:numel(simData(ii).entryTimesAfterAds)
            flagTime = (simData(ii).biolTimeVector >= simData(ii).entryTimesAfterAds(i_entry));
            simData(ii).insideVectorBiolTime(flagTime) = i_entry;
        end

        %With respect to the time of microfluidic perfusion
        %Entry time since microfluidic perfusion
        simData(ii).injPhageTime = simData(ii).entryTimesAfterAds + tAds;

        %Intracellular phage number, with respect to time since microfluidic perfusion
        simData(ii).flowTimeVector = biolTimeVector;
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

%% 1.2 SIMULATION OF CELL FATE (VANILLA AS IN PNAS 2021) - Building on top of the entry dynamics simulation
clc
close all

%Simulation setup
fateTimeVector = 0:0.1:65; %minutes
MOItosimulate = nVector;

%For the random number generator for the single-cell volume
convFac = (1e9*1e15)/(6.022e23); %Converting from number of molecules in um^3 to nanomolars
V0_mean = 1; %Mean of the cell length distribution
ss = 0;

%Threshold for decision making
CINorm = 56.4071;
CroNorm = 52.9038;
CIINorm = 221.9365;
lyticThresh = 21.59;
lysogenicThresh = 3.198;

%Parameters from Yao et al., PNAS 2021
fateFolder = [pwd '/Dependent files/'];
fateFileName = {[fateFolder 'exp2_fits.txt']};

ncol = 41;
bigNum = 1e5;
fateParams = zeros(bigNum, ncol);
for i_file = 1:length(fateFileName)
    fileID = fopen(fateFileName{i_file}, 'r');
    formatSpec = '%f'; %floating point numbers
    sizeA = [ncol Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    fateParams( (1 + (i_file-1)*size(A, 2)):(i_file*size(A, 2)), :) = A.';
    fclose(fileID);
end

%Delete zero rows, keep only solutions
failNum = 1e10;
fateParams(~any(fateParams,2), :) = [];
fateParams = fateParams(~all(fateParams == failNum, 2), :);
fateParams = sortrows(fateParams, ncol); %Sort solutions by last column

% Column 1-9: Production rates
%       1-3: Basal cI translation rate from PRM, PRM active/basal ratio, PRE/PRM basal activity ratio
%       4-6: cI translation rate, cro transcription rate, cro translation rate
%       7-9: cII transcription rate, cII translation rate, viral replication rate
%
% Column 10-16: Degradation rates
%       10-12: Dilution rate, cI mRNA degradation rate, CI degradation rate
%       13-14: cro mRNA degradation rate, Cro degradation rate
%       15-16: cII mRNA degradation rate, CII degradation rate
%
% Column 17-27: Hill coefficients
%       17-20: nPRM_CI (active), nPRM_CI (repress), nPRM_Cro, nPRE
%       21-27: ncro_Cro, ncro_CI, ncII_Cro, ncII_CI, nM_Cro, nM_CI, nDeg_CII
%
% Column 28-38: Activation/repression thresholds
%       28-31: KPRM_CI (active), KPRM_CI (repress), KPRM_Cro, KPRE
%       32-38: Kcro_Cro, Kcro_CI, KcII_Cro, KcII_CI, KM_Cro, KM_CI, KDeg_CII
%
% Column 39: Time offset between smFISH experiments and P+ qPCR experiment
% Column 40: Time offset for onset of replication
%
% Column 41: An objective function to evaluate the fitting

disp('Preparation for the cell-fate simulation is complete.');

%Options to solve ODEs
i_param = 1; %The principal parameter set only
options = odeset('Nonnegative', [], 'RelTol', 1e-6, 'AbsTol', 1e-6);
species = struct(); %Each row is a cell, all grouped within the same "row" of the master 'sim' data structure
outcomeEndtime = struct();

%Assign parameters
%cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
%CI translation rate, cro prod rate, Cro translation rate, cII prod rate
%CII translation rate, replication rate
prod = fateParams(i_param, 1:9);
prod(3) = prod(1)*prod(3);

%kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
degr = [fateParams(i_param, 10:16) 0];

%nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, nM_Cro, nM_CI, nDeg_CII
n = fateParams(i_param, 17:27);

%KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI,
%KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
K = fateParams(i_param, 28:38);
dt = fateParams(i_param, end-2); %offset
tau = fateParams(i_param, end-1);
kdil = degr(1);

%Simulation for each outside MOI
timerStart = tic;
fprintf('(Vanilla version) Now simulating for outside MOI =');
for i_moi = 1:numel(MOItosimulate)
    currentOutside = MOItosimulate(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    fprintf(' %d', currentOutside);

    for i_cell = idCellsWithMOI
        %Initial conditions for this particular cell, in concentration terms
        V0 = V0_mean;
        currentMOI = currentOutside;
        y0_ = [zeros(1, 6) currentMOI * convFac / V0];
        species.Vseries = V0.*exp(fateTimeVector.*kdil)';

        %-------------------- P- ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII]
        [tm_OP, ym_OP] = ode15s(@fv19, fateTimeVector, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), currentMOI, V0, convFac, ss);
        species.PminusConc = [tm_OP, ym_OP];
        species.PminusNum = [tm_OP, ym_OP .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_OP(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_OP(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pminus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pminus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pminus = 2;
        else
            outcomeEndtime.Pminus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pminus = 1;
            else
                outcomeFirst.Pminus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pminus = 1;
        elseif any(flagCI)
            outcomeFirst.Pminus = 2;
        else
            outcomeFirst.Pminus = 0;
        end

        %-------------------- P+ ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII] [lambda]
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, fateTimeVector, y0_, options, n, prod, ...
            degr, K, tau, V0_mean, convFac);
        species.PplusConc = [tm_wt, ym_wt];
        species.PplusNum = [tm_wt, ym_wt .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_wt(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_wt(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pplus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pplus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pplus = 2;
        else
            outcomeEndtime.Pplus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pplus = 1;
            else
                outcomeFirst.Pplus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pplus = 1;
        elseif any(flagCI)
            outcomeFirst.Pplus = 2;
        else
            outcomeFirst.Pplus = 0;
        end

        %Saving the results of this current cell to the master data structure
        simData(i_cell).totalInputMOI = currentOutside;
        simData(i_cell).fateTimeVector = fateTimeVector;
        simData(i_cell).defaultSize = V0;

        simData(i_cell).vanillaFate.species = species;
        simData(i_cell).vanillaFate.outcomeEndtime = outcomeEndtime;
        simData(i_cell).vanillaFate.outcomeFirst = outcomeFirst;
    end
end

timerEnd = toc(timerStart);

disp([newline 'Simulation and calculation are complete, taking ' num2str(timerEnd, '%.2f') ' sec.']);

%% 1.3 SIMULATION FOR CELL FATE (WITH LOGNORMAL-DISTRIBUTED CELL SIZE ONLY) -- Run the above first
clc
close all

%Additional setups
mu = V0_mean;
sigma = 0.28;
% mu = log(V0_mean) - sigma^2/2;
% volFunc = @(x) 1./(x .* sigma .* sqrt(2.*pi)) .* exp(-(log(x) - mu).^2./(2.*sigma.^2)); %PDF of log-normal distribution with mean = 1

%Simulation of cell size first for all cells in the dataset
for i_moi = 1:numel(MOItosimulate)
    currentOutside = MOItosimulate(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);
%     V0array = -slicesample(1, numel(idCellsWithMOI), 'pdf', volFunc, 'burnin', 1000, 'thin', 10, 'width', 1);
    V0array = lognrnd(log(mu), sigma, 1, numel(idCellsWithMOI));
    ii = 1;

    %Allocating random cell size to the sample
    for i_cell = idCellsWithMOI
        simData(i_cell).lognormalSize = V0array(ii);
        ii = ii+1;
    end
end

%Simulation for each outside MOI
timerStart = tic;
fprintf('(With log-normal cell size) Now simulating for outside MOI =');
for i_moi = 1:numel(MOItosimulate)
    currentOutside = MOItosimulate(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    fprintf(' %d', currentOutside);
    for i_cell = idCellsWithMOI
        %Initial conditions for this particular cell, in concentration terms
        V0 = simData(i_cell).lognormalSize;
        currentMOI = currentOutside;
        y0_ = [zeros(1, 6) currentMOI * convFac / V0];
        species.Vseries = V0.*exp(fateTimeVector.*kdil)';

        %-------------------- P- ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII]
        [tm_OP, ym_OP] = ode15s(@fv19, fateTimeVector, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), currentMOI, V0, convFac, ss);
        species.PminusConc = [tm_OP, ym_OP];
        species.PminusNum = [tm_OP, ym_OP .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_OP(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_OP(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pminus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pminus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pminus = 2;
        else
            outcomeEndtime.Pminus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pminus = 1;
            else
                outcomeFirst.Pminus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pminus = 1;
        elseif any(flagCI)
            outcomeFirst.Pminus = 2;
        else
            outcomeFirst.Pminus = 0;
        end

        %-------------------- P+ ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII] [lambda]
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, fateTimeVector, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        species.PplusConc = [tm_wt, ym_wt];
        species.PplusNum = [tm_wt, ym_wt .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_wt(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_wt(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pplus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pplus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pplus = 2;
        else
            outcomeEndtime.Pplus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pplus = 1;
            else
                outcomeFirst.Pplus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pplus = 1;
        elseif any(flagCI)
            outcomeFirst.Pplus = 2;
        else
            outcomeFirst.Pplus = 0;
        end

        %Saving the results of this current cell to the master data structure
        simData(i_cell).wsizeFate.species = species;
        simData(i_cell).wsizeFate.outcomeEndtime = outcomeEndtime;
        simData(i_cell).wsizeFate.outcomeFirst = outcomeFirst;
    end
end

timerEnd = toc(timerStart);

disp([newline 'Simulation and calculation are complete, taking ' num2str(timerEnd, '%.2f') ' sec.']);

%% 1.4 SIMULATION FOR CELL FATE (WITH SIMULATED ENTRY DYNAMICS ONLY) -- Run the above first
clc
close all

%Additional setups
tDeadline = 5; %min, the deadline after which further phage entries do not contribute to the decision

%Simulation for each outside MOI
timerStart = tic;
fprintf('(With entry dynamics, by 5-min deadline) Now simulating for outside MOI =');
for i_moi = 1:numel(MOItosimulate)
    currentOutside = MOItosimulate(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    fprintf(' %d', currentOutside);
    for i_cell = idCellsWithMOI
        %Initial conditions for this particular cell, in concentration terms
        V0 = V0_mean;
        currentMOI = sum([simData(i_cell).entryTimesAfterAds] <= (tDeadline*60));
        y0_ = [zeros(1, 6) currentMOI * convFac / V0];
        species.Vseries = V0.*exp(fateTimeVector.*kdil)';

        %-------------------- P- ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII]
        [tm_OP, ym_OP] = ode15s(@fv19, fateTimeVector, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), currentMOI, V0, convFac, ss);
        species.PminusConc = [tm_OP, ym_OP];
        species.PminusNum = [tm_OP, ym_OP .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_OP(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_OP(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pminus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pminus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pminus = 2;
        else
            outcomeEndtime.Pminus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pminus = 1;
            else
                outcomeFirst.Pminus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pminus = 1;
        elseif any(flagCI)
            outcomeFirst.Pminus = 2;
        else
            outcomeFirst.Pminus = 0;
        end

        %-------------------- P+ ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII] [lambda]
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, fateTimeVector, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        species.PplusConc = [tm_wt, ym_wt];
        species.PplusNum = [tm_wt, ym_wt .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_wt(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_wt(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pplus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pplus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pplus = 2;
        else
            outcomeEndtime.Pplus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pplus = 1;
            else
                outcomeFirst.Pplus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pplus = 1;
        elseif any(flagCI)
            outcomeFirst.Pplus = 2;
        else
            outcomeFirst.Pplus = 0;
        end

        %Saving the results of this current cell to the master data structure
        simData(i_cell).effectiveInputMOI = currentMOI;
        
        simData(i_cell).wentryFate.species = species;
        simData(i_cell).wentryFate.outcomeEndtime = outcomeEndtime;
        simData(i_cell).wentryFate.outcomeFirst = outcomeFirst;
    end
end

timerEnd = toc(timerStart);

disp([newline 'Simulation and calculation are complete, taking ' num2str(timerEnd, '%.2f') ' sec.']);

%% 1.5 SIMULATION FOR CELL FATE (WITH BOTH CELL SIZE AND ENTRY DYNAMICS) -- Run the above first
clc
close all

%No additional setups

%Simulation for each outside MOI
timerStart = tic;
fprintf('(With log-normal cell size and 5-min entry deadline) Now simulating for outside MOI =');
for i_moi = 1:numel(MOItosimulate)
    currentOutside = MOItosimulate(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    fprintf(' %d', currentOutside);
    for i_cell = idCellsWithMOI
        %Initial conditions for this particular cell, in concentration terms
        V0 = simData(i_cell).lognormalSize;
        currentMOI = numel(simData(i_cell).entryTimesAfterAds <= (tDeadline*60));
        y0_ = [zeros(1, 6) currentMOI * convFac / V0];
        species.Vseries = V0.*exp(fateTimeVector.*kdil)';

        %-------------------- P- ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII]
        [tm_OP, ym_OP] = ode15s(@fv19, fateTimeVector, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), currentMOI, V0, convFac, ss);
        species.PminusConc = [tm_OP, ym_OP];
        species.PminusNum = [tm_OP, ym_OP .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_OP(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_OP(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pminus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pminus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pminus = 2;
        else
            outcomeEndtime.Pminus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pminus = 1;
            else
                outcomeFirst.Pminus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pminus = 1;
        elseif any(flagCI)
            outcomeFirst.Pminus = 2;
        else
            outcomeFirst.Pminus = 0;
        end

        %-------------------- P+ ------------------------------------------
        %Columns of output: [cI mRNA] [cro mRNA] [cII mRNA] [CI] [Cro] [CII] [lambda]
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, fateTimeVector, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        species.PplusConc = [tm_wt, ym_wt];
        species.PplusNum = [tm_wt, ym_wt .* species.Vseries ./ convFac];

        %Decision-making
        flagCro = (ym_wt(:, 5)/CroNorm) > lyticThresh;
        flagCI = (ym_wt(:, 4)/CINorm) > lysogenicThresh;

        %End-time decision
        if any(flagCro) && any(flagCI)
            outcomeEndtime.Pplus = 3;
        elseif any(flagCro)
            outcomeEndtime.Pplus = 1;
        elseif any(flagCI)
            outcomeEndtime.Pplus = 2;
        else
            outcomeEndtime.Pplus = 0;
        end

        %Fist-passage decision
        if any(flagCro) && any(flagCI)
            if find(flagCro, 1) <= find(flagCI, 1)
                outcomeFirst.Pplus = 1;
            else
                outcomeFirst.Pplus = 2;
            end
        elseif any(flagCro)
            outcomeFirst.Pplus = 1;
        elseif any(flagCI)
            outcomeFirst.Pplus = 2;
        else
            outcomeFirst.Pplus = 0;
        end

        %Saving the results of this current cell to the master data structure
        simData(i_cell).fullFate.species = species;
        simData(i_cell).fullFate.outcomeEndtime = outcomeEndtime;
        simData(i_cell).fullFate.outcomeFirst = outcomeFirst;
    end
end

timerEnd = toc(timerStart);

disp([newline 'Simulation and calculation are complete, taking ' num2str(timerEnd, '%.2f') ' sec.']);

%% 2. PREPARATIONS FOR VISUALIZATION
clc
close all

%Loading Lanying's data from Cell 2010, from Figure 2C
LZ.moi = [0 1 2 3 4 5];
LZ.mean = [0 33.95248568 55.20241333 64.96599919 64.59430523 73.68221541]/100;
LZ.error = [0 35.84423288 58.5809364 70.23626946 70.67538631 82.33086405]/100;
LZ.color = [70 140 210]/255;

%Colors of MOI as in Cell 2010, with blank for MOI = 0 for completeness
moiColors = [0 0 0; 40 30 125; 185 30 40; 45 150 70; 30 150 230; 185 5 120]/255;

%Fate colors and names
lysogenicColor = [250 80 80]/255;
lyticColor = [70 200 90]/255;
fatesColor = [0.5 0.5 0.5; lyticColor; lysogenicColor; 0.8 0.7 0.1];
fateNames = ["Failed" "Lytic" "Lysogenic" "Mixed"]; %0, 1, 2, and 3

%For histograms later
ticksMOIVector = 0:1:max(nVector);
binningMOIVector = (-0.5):1:(max(nVector)+0.5);

%From the simData structure
modelNames = ["Vanilla" "With log-normal size" "With entry dynamics" "With size and entry"];
modelFieldnames = ["vanillaFate" "wsizeFate" "wentryFate" "fullFate"];
modelSizeFields = ["defaultSize" "lognormalSize" "defaultSize" "lognormalSize"];
modelMOIFields = ["totalInputMOI" "totalInputMOI" "effectiveInputMOI" "effectiveInputMOI"];

%Defining the model
mSpace = 0:0.05:10;
Hill.parameters = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Lower', [1 0 0 0], ...
    'Upper', [1 0 50 5], ...
    'Startpoint', [1 0 1 1]);
Hill.function = @(a, b, h, K, m) a .* m.^h ./ (m.^h + K.^h) + b;
Hill.model = fittype(Hill.function, 'Options', Hill.parameters,...
    'Independent', 'm', 'Coefficient', {'a', 'b', 'h', 'K'});

disp('Preparation is complete.');

%% 3.1 DISTRIBUTION OF CELL SIZE
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);

%Theoretical calculation
pdf_lognormal = @(x) 1./(x .* sigma .* sqrt(2.*pi)) .* exp(-(log(x) - log(mu)).^2./(2.*sigma.^2));
xspace = 0:0.01:4;
yspace = pdf_lognormal(xspace);

for i_moi = 1:numel(nVector)
    subplot(2, 3, i_moi);
    hold on;

    currentOutside = nVector(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    graphSim = histogram([simData(idCellsWithMOI).lognormalSize], 0:0.1:3,...
        'Normalization', 'PDF', 'FaceColor', [1 1 1]*0.6,...
        'DisplayName', ['Simulation (\itN\rm = ' num2str(numCellsPerMOI) ' cells)']);

    graphAna = plot(xspace, yspace, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5,...
        'DisplayName', ['Log-normal PDF: \mu = ' num2str(mu) ', \sigma = ' num2str(sigma)]);

    %Adjusting plot properties
    legend([graphSim, graphAna], 'Location', 'NorthEast', 'FontSize', 10)
    xlim([0 2.5]);
    xlabel('Cell length (\mum)   |   Cell volume (\mum^3)');
    ylabel('Probability density');
    ylim([0 max(get(gca, 'YLim'))*1.4]);
    title(['For outside MOI = ' num2str(currentOutside)], 'FontSize', 14);
    grid on; box on;
end

%% 3.2 DISTRIBUTION OF DIFFERENT INPUT MOI
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);

for i_moi = 1:numel(nVector)
    currentOutside = nVector(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

    subplot(3, numel(nVector), i_moi);
    hold on;

    histogram([simData(idCellsWithMOI).totalInputMOI], binningMOIVector,...
        'Normalization', 'Probability', 'FaceColor', [1 1 1]*0.6);

    %Adjusting plot properties
    xticks(ticksMOIVector);
    xlabel('Inside MOI');
    ylabel('Frequency');
    ylim([0 max(get(gca, 'YLim'))*1.1]);
    title(['Outside MOI = ' num2str(currentOutside) ', if 100% entry'], 'FontSize', 10);
    grid on; box on;

    subplot(3, numel(nVector), numel(nVector) + i_moi);
    hold on;

    histogram([simData(idCellsWithMOI).insideMaxCount], binningMOIVector,...
        'Normalization', 'Probability', 'FaceColor', [1 1 1]*0.6);

    %Adjusting plot properties
    xticks(ticksMOIVector);
    xlabel('Inside MOI');
    ylabel('Frequency');
    ylim([0 max(get(gca, 'YLim'))*1.1]);
    title(['For outside MOI = ' num2str(currentOutside) ', at \infty'], 'FontSize', 10);
    grid on; box on;

    subplot(3, numel(nVector), numel(nVector)*2 + i_moi);
    hold on;

    histogram([simData(idCellsWithMOI).effectiveInputMOI], binningMOIVector,...
        'Normalization', 'Probability', 'FaceColor', [1 1 1]*0.6);

    %Adjusting plot properties
    xticks(ticksMOIVector);
    xlabel('Inside MOI');
    ylabel('Frequency');
    ylim([0 max(get(gca, 'YLim'))*1.1]);
    title(['For outside MOI = ' num2str(currentOutside) ' at 5 min'], 'FontSize', 10);
    grid on; box on;
end

%% 4.1 CI AND CRO TRAJECTORIES IN DIFFERENT MODEL VERSIONS
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05+0 0.01 0.9 0.98]);

numExamples = 50;

for i_model = 1:numel(modelNames)
    for i_moi = 1:numel(nVector)
    currentOutside = nVector(i_moi);
    idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);
    idCellsToShow = idCellsWithMOI(randsample(numel(idCellsWithMOI), numExamples, 'false')); %Random sampling without replacements

    subplot(numel(modelNames), numel(nVector), (i_model-1)*numel(nVector) + i_moi);
    hold on;

    %Example trajectories
    for i_cell = idCellsToShow
        currentCell = simData(i_cell).(modelFieldnames(i_model));

        CInormalized = currentCell.species.PplusConc(:, 5) / CINorm;
        Cronormalized = currentCell.species.PplusConc(:, 6) / CroNorm;

        graphCell = plot(CInormalized, Cronormalized,...
            'LineStyle', '-', 'Color', [[1 1 1]*0 1], 'LineWidth', 1);
    end

    %Decision threshold
    xline(lysogenicThresh, 'Color', lysogenicColor, 'LineStyle', '-', 'LineWidth', 1.5);
    yline(lyticThresh, 'Color', lyticColor, 'LineStyle', '-', 'LineWidth', 1.5);
    
    %Adjusting plot properties
    axis square;
    axis([0 15 0 30]);
    xlabel('[CI]');
    ylabel('[Cro]');
    if i_model == 1
        title(['Outside MOI = ' num2str(currentOutside) newline,...
            char(modelNames(i_model))], 'FontSize', 10);
    else
        title([char(modelNames(i_model))], 'FontSize', 10);
    end
    grid on; box on;
    end
end

%% 4.2 DISTRIBUTION OF INFECTION OUTCOMES
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05+0 0.05 0.9 0.9]);
sgtitle(['Infection outcome by end time (65 min)']);

for i_model = 1:numel(modelNames)
    fateStats = zeros(numel(nVector), 4);
    for i_moi = 1:numel(nVector)
        currentOutside = nVector(i_moi);
        idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

        %Data compilation and statistics calculation
        vectorFate = [];
        for i_cell = idCellsWithMOI
            vectorFate = [vectorFate simData(i_cell).(modelFieldnames(i_model)).outcomeEndtime.Pplus];
        end
        for i_fate = [0 1 2 3] %Failed | Lytic | Lysogenic | Mixed
            fateStats(i_moi, i_fate+1) = sum(vectorFate == i_fate) / numel(vectorFate);
        end

        %Visualization
        subplot(numel(modelNames), numel(nVector), (i_model-1)*numel(nVector) + i_moi);
        hold on;

        graphFates = bar(1:4, fateStats(i_moi, :), 'FaceColor', 'Flat');

        %Adjusting plot properties
        xlim([0 5]);
        xticks(1:4);
        xticklabels(fateNames);
        ylim([0 max(get(gca, 'YLim'))*1.1]);
        ylabel('Frequency');
        if i_model == 1
            title(['Outside MOI = ' num2str(currentOutside) newline,...
                char(modelNames(i_model))], 'FontSize', 10);
        else
            title([char(modelNames(i_model))], 'FontSize', 10);
        end
        for i_fate = 1:4 %Recoloring the bars
            graphFates.CData(i_fate,:) = fatesColor(i_fate, :);
        end
        grid on; box on;
    end
end

%% 4.3 MOI RESPONSE CURVES (QUICKL HILL FITTING)
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05 0.25 0.9 0.5]);
sgtitle(['Between lytic and lysogenic outcomes, by end time (65 min)']);

for i_model = 1:numel(modelNames)
    %Fraction of lysogeny vs. MOI
    subplot(1, numel(modelNames), i_model);
    hold on;

    %Data compilation and statistics calculation
    fractionLysogeny = zeros(1, numel(nVector));
    vectorFate = cell(1, numel(nVector));
    for i_moi = 1:numel(nVector)
        currentOutside = nVector(i_moi);
        idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

        vectorFate{i_moi} = [];
        for i_cell = idCellsWithMOI
            vectorFate{i_moi} = [vectorFate{i_moi} simData(i_cell).(modelFieldnames(i_model)).outcomeEndtime.Pplus];
        end
        fractionLysogeny(i_moi) = sum(vectorFate{i_moi} == 2) / sum(vectorFate{i_moi} == 1 | vectorFate{i_moi} == 2);
    end

    %Fitting and graphing
    [xToFit, yToFit] = prepareCurveData(nVector, fractionLysogeny);
    [param, gof] = fit(xToFit, yToFit, Hill.model);
    fSpaceModel = Hill.function(param.a, param.b, param.h, param.K, mSpace);
    plot(mSpace, fSpaceModel, 'LineStyle', '-', 'Color', lysogenicColor, 'LineWidth', 1.5);

    graphModel = plot(nVector, fractionLysogeny,...
        'Marker', 'o', 'MarkerFaceColor', lysogenicColor, 'MarkerEdgeColor', 'None', 'LineStyle', 'None',...
        'DisplayName', ['\bfModel\rm: Hill = ' num2str(param.h, '%.2f') ', midpoint = ' num2str(param.K, '%.2f')]);

    [xToFit, yToFit] = prepareCurveData(LZ.moi, LZ.mean);
    [param, gof] = fit(xToFit, yToFit, Hill.model);
    fSpaceZeng = Hill.function(param.a, param.b, param.h, param.K, mSpace);
    plot(mSpace, fSpaceZeng, 'LineStyle', '-', 'Color', LZ.color, 'LineWidth', 1.5);

    graphLanying = errorbar(LZ.moi, LZ.mean, LZ.error-LZ.mean,...
        'Marker', 'v', 'MarkerFaceColor', LZ.color, 'MarkerEdgeColor', 'None',...
        'LineStyle', 'None', 'Color', LZ.color,...
        'DisplayName', ['\bfCell 2010\rm: Hill = ' num2str(param.h, '%.2f') ', midpoint = ' num2str(param.K, '%.2f')]);

    %Adjusting plot properties
    legend([graphModel graphLanying], 'Location', 'SouthOutside', 'FontSize', 10);
    xlim([0 5.5]);
    xlabel('Outside MOI, \itn');
    ylim([0 1.1]);
    ylabel('Fraction, \itf');
    title([modelNames(i_model)], 'FontSize', 12);
    grid on; box on;
end

%% 4.4 MOI RESPONSE CURVES (HILL FIT WITH BOOTSTRAPPING)
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.05 0.05 0.9 0.9]);
sgtitle(['Between lytic and lysogenic outcomes, by end time (65 min)']);

nHillParams = 4; %a, b, h, K
nSampling = 1000;

%Lanying's data first
subplot(2, 3, 1);
hold on;

warning('off');
%Bootstrapped fitting on the mean level
bootstrapParams = zeros(2, nHillParams);
paramSampling = zeros(nSampling, nHillParams);
for i_sampling = 1:nSampling
    [~, idx] = datasample(LZ.moi, numel(LZ.moi));
    [xToFit, yToFit] = prepareCurveData(LZ.moi(idx), LZ.mean(idx));
    [param, gof] = fit(xToFit, yToFit, Hill.model);
    paramSampling(i_sampling, :) = [param.a param.b param.h param.K];
end
bootstrapParams(1, :) = mean(paramSampling, 1);
bootstrapParams(2, :) = std(paramSampling, 1);
warning('on');

fSpaceZeng = Hill.function(bootstrapParams(1, 1), bootstrapParams(1, 2),...
    bootstrapParams(1, 3), bootstrapParams(1, 4), mSpace);
plot(mSpace, fSpaceZeng, 'LineStyle', '-', 'Color', LZ.color, 'LineWidth', 1.5);

graphLanying = errorbar(LZ.moi, LZ.mean, LZ.error-LZ.mean,...
    'Marker', 'v', 'MarkerFaceColor', LZ.color, 'MarkerEdgeColor', 'None',...
    'LineStyle', 'None', 'Color', LZ.color,...
    'DisplayName', ['Hill coeff = ' num2str(bootstrapParams(1, 3), '%.2f') ' \pm ' num2str(bootstrapParams(2, 3), '%.2f') newline,...
    'Midpoint = ' num2str(bootstrapParams(1, 4), '%.2f') ' \pm ' num2str(bootstrapParams(2, 4), '%.2f')]);

%Adjusting plot properties
legend([graphLanying], 'Location', 'SouthEast', 'FontSize', 10);
xlim([0 5.5]);
xlabel('Outside MOI, \itn');
ylim([0 1.1]);
ylabel('Fraction, \itf');
title(['Zeng et al., Cell 2010'], 'FontSize', 12);
grid on; box on;

%Saving for later
hillCoeff(1, :) = [bootstrapParams(1, 3) bootstrapParams(2, 3)];

%For TN's model predictions
for i_model = 1:numel(modelNames)
    %Fraction of lysogeny vs. MOI
    subplot(2, 3, i_model+1);
    hold on;

    %Data compilation and statistics calculation
    fractionLysogeny = zeros(1, numel(nVector));
    vectorFate = cell(1, numel(nVector));
    for i_moi = 1:numel(nVector)
        currentOutside = nVector(i_moi);
        idCellsWithMOI = find([simData(:).outsideMaxCount] == currentOutside);

        vectorFate{i_moi} = [];
        for i_cell = idCellsWithMOI
            vectorFate{i_moi} = [vectorFate{i_moi} simData(i_cell).(modelFieldnames(i_model)).outcomeEndtime.Pplus];
        end
        fractionLysogeny(i_moi) = sum(vectorFate{i_moi} == 2) / sum(vectorFate{i_moi} == 1 | vectorFate{i_moi} == 2);
    end

    warning('off');
    %Bootstrapped fitting on the single-cell level
    bootstrapParams = zeros(2, nHillParams);
    paramSampling = zeros(nSampling, nHillParams);
    for i_sampling = 1:nSampling
        [~, idx] = datasample(vectorFate{i_moi}, numel(vectorFate{i_moi}));
        bootstrapFraction = zeros(1, numel(nVector));
        for i_moi = 1:numel(nVector)
            bootstrapFraction(i_moi) = sum(vectorFate{i_moi}(idx) == 2) / sum(vectorFate{i_moi}(idx) == 1 | vectorFate{i_moi}(idx) == 2);
        end
        [xToFit, yToFit] = prepareCurveData(nVector, bootstrapFraction);
        [param, gof] = fit(xToFit, yToFit, Hill.model);
        paramSampling(i_sampling, :) = [param.a param.b param.h param.K];
    end
    bootstrapParams(1, :) = mean(paramSampling, 1);
    bootstrapParams(2, :) = std(paramSampling, 1);
    warning('on');

    fSpace = Hill.function(bootstrapParams(1, 1), bootstrapParams(1, 2),...
        bootstrapParams(1, 3), bootstrapParams(1, 4), mSpace);
    plot(mSpace, fSpace, 'LineStyle', '-', 'Color', lysogenicColor, 'LineWidth', 1.5);

    graphModel = plot(nVector, fractionLysogeny,...
        'Marker', 'o', 'MarkerFaceColor', lysogenicColor, 'MarkerEdgeColor', 'None', 'LineStyle', 'None',...
        'DisplayName', ['Hill coeff = ' num2str(bootstrapParams(1, 3), '%.2f') ' \pm ' num2str(bootstrapParams(2, 3), '%.2f') newline,...
        'Midpoint = ' num2str(bootstrapParams(1, 4), '%.2f') ' \pm ' num2str(bootstrapParams(2, 4), '%.2f')]);

    %Adjusting plot properties
    legend([graphModel], 'Location', 'SouthEast', 'FontSize', 10);
    xlim([0 5.5]);
    xlabel('Outside MOI, \itn');
    ylim([0 1.1]);
    ylabel('Fraction, \itf');
    title([modelNames(i_model)], 'FontSize', 12);
    grid on; box on;

    %Saving for later
    hillCoeff(i_model+1, :) = [bootstrapParams(1, 3) bootstrapParams(2, 3)];
end

%% 4.5 COMPARING HILL COEFFICIENTS
clc
close all

figure('Unit', 'Normalized', 'OuterPosition', [0.25 0.05 0.5 0.9]);
hold on;

modelsInterest = [3 5 1];
bar(1:3, hillCoeff(modelsInterest, 1), 'FaceColor', lysogenicColor);
errorbar(1:3, hillCoeff(modelsInterest, 1), hillCoeff(modelsInterest, 2),...
    'Marker', 'None', 'LineStyle', 'None', 'LineWidth', 1.2, 'CapSize', 12, 'Color', 'k');

%Adjusting plot properties
axis square;
xlim([0 4]);
xticks(1:3);
xticklabels([["PNAS 2021\newline+ Variable size"],...
    ["PNAS 2021\newline+ Variable size\newline+ Entry dynamics"],...
    ["Data from Cell 2010\newline(refitted)"]]);
ylabel(['Coefficient \pm SE by bootstrapping']);
title(['Hill coefficients of the MOI response curves'], 'FontSize', 12);
grid on; box on;

