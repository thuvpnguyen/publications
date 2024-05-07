%% Overall comments:
%This code include the following parts:
%1. Seth load fits: read the fitted data from exp2_fits.txt. Must run, unless refit (which I cannot now)
%2. Parameter statistics (5 sections; analysis that does not depend on following simulations)
%3. Simulations (3 sections, including Seth's simulation and my simulation)
%4. Figure creation (including Seth's Fig3 and my various trials)


%% Seth load fits
clear
clc

filenames = {'exp2_fits.txt'};

ncol = 41;
bigNum = 1e5;
data = zeros(bigNum, ncol);
for i = 1:length(filenames)
    fileID = fopen(filenames{i}, 'r');
    formatSpec = '%f'; %floating point numbers
    sizeA = [ncol Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    data( (1 + (i-1)*size(A, 2)):(i*size(A, 2)), :) = A.';
    fclose(fileID);
end

%Delete zero rows, keep only solutions
failNum = 1e10;
data(~any(data,2), :) = [];
data = data(~all(data == failNum, 2), :);

%Sort solutions by last column
data = sortrows(data, ncol);

%% TY comments: every column refers to a parameter
%Column 1-9: Production rates
%       1-3: Basal cI translation rate from PRM, PRM active/basal ratio, PRE/PRM basal activity ratio
%       4-6: cI translation rate, cro transcription rate, cro translation rate
%       7-9: cII transcription rate, cII translation rate, viral replication rate

%Column 10-16: Degradation rates
%       10-12: Dilution rate, cI mRNA degradation rate, CI degradation rate
%       13-14: cro mRNA degradation rate, Cro degradation rate
%       15-16: cII mRNA degradation rate, CII degradation rate

%Column 17-27: Hill coefficients
%       17-20: nPRM_CI (active), nPRM_CI (repress), nPRM_Cro, nPRE
%       21-27: ncro_Cro, ncro_CI, ncII_Cro, ncII_CI, nM_Cro, nM_CI, nDeg_CII
    
%Column 28-38: Activation/repression thresholds
%       28-31: KPRM_CI (active), KPRM_CI (repress), KPRM_Cro, KPRE
%       32-38: Kcro_Cro, Kcro_CI, KcII_Cro, KcII_CI, KM_Cro, KM_CI, KDeg_CII

%Column 39: Time offset between smFISH experiments and P+ qPCR experiment
%Column 40: Time offset for onset of replication

%Column 41: An objective function to evaluate the fitting

%% TY parameter statistics (for Table S7)
text = {'K_{tx}^{1, RM}', '\alpha_{RM}', '\alpha_{RE}', 'k_{tr}^1', ...    1-4
        'k_{tx}^2', 'k_{tr}^2', 'k_{tx}^3', 'k_{tr}^3', 'k_{\lambda}', ... 5-9
        'k_d', 'k_m^1', 'k_P^1', 'k_m^2', 'k_P^2', 'k_m^3', 'k_P^3', ...   10-16
        'n_{RM,CI}^a', 'n_{RM,CI}^r', 'n_{RM,Cro}', 'n_{RE}', ...          17-20
        'n_{cro,Cro}', 'n_{cro,CI}', 'n_{cII,Cro}', 'n_{cII,CI}', ...      21-24
        'n_{\lambda,Cro}', 'n_{\lambda,CI}', 'n_{CII}', ...                25-27
        'K_{RM,CI}^a', 'K_{RM,CI}^r', 'K_{RM,Cro}', 'K_{RE}', ...          28-31
        'K_{cro,Cro}', 'K_{cro,CI}', 'K_{cII,Cro}', 'K_{cII,CI}', ...      32-35
        'K_{\lambda,Cro}', 'K_{\lambda,CI}', 'K_{CII}', ...                36-38
        '\deltat', '\tau_{\lambda}'};                                     %39-40

row2c = 1:44;       %Only look at the best 44 fits

%Cook the data to perform all the normalizations
data_ = zeros(max(row2c), 40);
desc_ = cell(40, 1);

for ipar = 1:40
    if ismember(ipar, [4 28 29 33 35 37])    %normalized by Kcro_CI (33)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 33);
        desc_{ipar} = [text{ipar} '/' text{33}];
    elseif ismember(ipar, [6 30 32 34 36])   %normalized by KPRM_Cro (30)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 30);
        desc_{ipar} = [text{ipar} '/' text{30}];
    elseif ismember(ipar, [8 31 38])         %normalized by KPRE (31)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 31);
        desc_{ipar} = [text{ipar} '/' text{31}];
    else                                     %no normalization
        data_(:, ipar) = data(row2c, ipar);
        desc_{ipar} = text{ipar};
    end
end

%Table S7
output_ = zeros(40, 5);
for iparam = 1:40
    %Ensemble fit value (main text)
    output_(iparam, 1) = round(data_(1, iparam), 2, 'significant');
    
    %Ensemble mean and standard deviation
    output_(iparam, 2) = round(mean(data_(:, iparam)), 2, 'significant');
    output_(iparam, 3) = round(std(data_(:, iparam)), 2, 'significant');
    
    %Ensemble range
    output_(iparam, 4) = round(min(data_(:, iparam)), 2, 'significant');
    output_(iparam, 5) = round(max(data_(:, iparam)), 2, 'significant');
end

%% TY test step 1: simulate all 44 fits 
row2c = 1:44;
% row2c = 1;
clear sol
sol = cell(numel(row2c), 1);

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
ss = 0;

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, 'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

for ii = 1:numel(row2c)
    sol{ii} = struct();
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod = data(ii, 1:9);
    prod(3) = prod(1)*prod(3);

    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    degr = [data(ii, 10:16) 0];
    
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, nM_Cro, nM_CI, nDeg_CII
    n = data(ii, 17:27);
    
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = data(ii, 28:38);
    
    dt = data(ii, end-2); %offset
    tau = data(ii, end-1);
    kdil = degr(1);
    sol{ii}.parameters = data(ii, :);
    sol{ii}.prod = prod;
    sol{ii}.degr = degr;
    sol{ii}.n = n;
    sol{ii}.K = K;
    sol{ii}.dt = dt;
    sol{ii}.tau = tau;
    
    V = V0.*exp(tspan.*kdil)';
    
    for imoi = 1:5
        y0_ = [zeros(1, 6) imoi * convFac / V0];
        %-------------------- P- ------------------------------------------
        name_1 = ['m' num2str(imoi)];
        name_2 = ['m' num2str(imoi) 'Num'];
        
        [tm_OP, ym_OP] = ode15s(@fv19, tspan, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), imoi, V0, convFac, ss);
        sol{ii}.(name_1) = [tm_OP, ym_OP];
        sol{ii}.(name_2) = [tm_OP, ym_OP .* V ./ convFac];
        
        %-------------------- P+ ------------------------------------------
        name_3 = ['m' num2str(imoi) '_wt'];
        name_4 = ['m' num2str(imoi) 'Num_wt'];
        
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, tspan, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        sol{ii}.(name_3) = [tm_wt, ym_wt];
        sol{ii}.(name_4) = [tm_wt, ym_wt .* V ./ convFac];
        
    end
end

%% TY test step 2: create CI-Cro trajectories (norm by max, norm by thresh)
cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;
% cm = [255   0   0;
%       101 185 222;
%        75 123 190;
%        88  78 160;
%         7  24  50]./255;
    
fig3_new = figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.5 0.9], ...
    'name', 'Fig 3: Phase plane traj.');    

for ii = 1:numel(row2c)
    %CI, Cro, CII normalizations
    CINorm = sol{ii}.K(6); %KCro_CI
    CroNorm = sol{ii}.K(3); %KPRM_Cro
    CIINorm = sol{ii}.K(4); %KPRE
    
    CI_Pm_lst = zeros(numel(tspan), 5);
    Cro_Pm_lst = zeros(numel(tspan), 5);
    CI_Pp_lst = zeros(numel(tspan), 5);
    Cro_Pp_lst = zeros(numel(tspan), 5);
    for imoi = 1:5
        data_Pm = sol{ii}.(['m' num2str(imoi)]);
        CI_Pm = data_Pm(:, 5)/CINorm;      %CI
        Cro_Pm = data_Pm(:, 6)/CroNorm;    %Cro
        CI_Pm_lst(:, imoi) = CI_Pm;
        Cro_Pm_lst(:, imoi) = Cro_Pm;
        
        data_Pp = sol{ii}.(['m' num2str(imoi) '_wt']);
        CI_Pp = data_Pp(:, 5)/CINorm;      %CI
        Cro_Pp = data_Pp(:, 6)/CroNorm;    %Cro
        CI_Pp_lst(:, imoi) = CI_Pp;
        Cro_Pp_lst(:, imoi) = Cro_Pp;
    end
    %find maximal for CI and Cro in every fit
    CI_Pm_mx = max(CI_Pm_lst(:, 5));
    Cro_Pm_mx = max(Cro_Pm_lst(:, 5));
    CI_Pp_mx = max(CI_Pp_lst(:, 5));
    Cro_Pp_mx = max(Cro_Pp_lst(:, 5));
    
    %define lytic and lysogenic thresholds in every fit
    sol{ii}.lytic_thresh_max = max(Cro_Pp_lst(:, 1));
    sol{ii}.lytic_thresh_min = max(Cro_Pm_lst(:, 5));
    sol{ii}.lysog_thresh_max = max(CI_Pp_lst(:, 2));
    sol{ii}.lysog_thresh_min = max(CI_Pp_lst(:, 1));
    lytic_perc = 0.5;        %0: lower bound; 1: higher bound
    lysog_perc = 0.5;
    
    lytic_thresh = sol{ii}.lytic_thresh_min * (1-lytic_perc) + sol{ii}.lytic_thresh_max * lytic_perc;
    lysog_thresh = sol{ii}.lysog_thresh_min * (1-lysog_perc) + sol{ii}.lysog_thresh_max * lysog_perc;
    
%     sol{ii}.lytic_thresh = mean([lytic_thresh_max lytic_thresh_min]);
%     sol{ii}.lysog_thresh = mean([lysog_thresh_max lysog_thresh_min]);
    
    for imoi = 1:5
        %Phase plane trajectory -------------------------------------------
        subplot(1, 2, 1); hold on;      %P-
        
%         %no normalization
%         xplot = CI_Pm_lst(:, imoi);
%         yplot = Cro_Pm_lst(:, imoi);
%         norm_txt = '';
        
%         %normalize with max CI & Cro
%         xplot = CI_Pm_lst(:, imoi) / CI_Pm_mx;
%         yplot = Cro_Pm_lst(:, imoi) / Cro_Pm_mx;
%         norm_txt = ' by max at MOI=5';
        
        %normalize with lytic and lysogenic thresholds
        xplot = CI_Pm_lst(:, imoi) / lysog_thresh;
        yplot = Cro_Pm_lst(:, imoi) / lytic_thresh;
        norm_txt = ' by decision threshold';
        
        p_ = plot(xplot, yplot, '-', 'color', cm(imoi, :), ...
            'linewidth', 0.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
        if imoi == 5
            xlabel(['CI concentration (normalized' norm_txt ')']);
            ylabel(['Cro concentration (normalized' norm_txt ')']);
            xlim([0 3]);
            ylim([0 1.5]);
%             pbaspect([1, 3, 1]);
%             pbaspect([1, 2, 1]);
            title('Cro/CI trajectory at P-');
            pbaspect([1, 1, 1]);
        end
        
        subplot(1, 2, 2); hold on;      %P+
        
%         %no normalization
%         xplot = CI_Pm_lst(:, imoi);
%         yplot = Cro_Pm_lst(:, imoi);
%         norm_txt = '';
        
%         %normalize with max CI & Cro
%         xplot = CI_Pp_lst(:, imoi) / CI_Pp_mx;
%         yplot = Cro_Pp_lst(:, imoi) / Cro_Pp_mx;
%         norm_txt = ' by max at MOI=5';
        
        %normalize with lytic and lysogenic thresholds
        xplot = CI_Pp_lst(:, imoi) / lysog_thresh;
        yplot = Cro_Pp_lst(:, imoi) / lytic_thresh;
        norm_txt = ' by decision threshold';
        
        p_ = plot(xplot, yplot, '-', 'color', cm(imoi, :), ...
            'linewidth', 0.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        if imoi == 5
            xlabel(['CI concentration (normalized' norm_txt ')']);
            ylabel(['Cro concentration (normalized' norm_txt ')']);
            xlim([0 3]);
            ylim([0 1.5]);
%             pbaspect([1, 3, 1]);
%             pbaspect([1, 2, 1]);
            title('Cro/CI trajectory at P+');
            pbaspect([1, 1, 1]);
        end
    end
end


%% TY parameter statistics (for Table S10)
text = {'K_{tx}^{1, RM}', '\alpha_{RM}', '\alpha_{RE}', 'k_{tr}^1', ...    1-4
        'k_{tx}^2', 'k_{tr}^2', 'k_{tx}^3', 'k_{tr}^3', 'k_{\lambda}', ... 5-9
        'k_d', 'k_m^1', 'k_P^1', 'k_m^2', 'k_P^2', 'k_m^3', 'k_P^3', ...   10-16
        'n_{RM,CI}^a', 'n_{RM,CI}^r', 'n_{RM,Cro}', 'n_{RE}', ...          17-20
        'n_{cro,Cro}', 'n_{cro,CI}', 'n_{cII,Cro}', 'n_{cII,CI}', ...      21-24
        'n_{\lambda,Cro}', 'n_{\lambda,CI}', 'n_{CII}', ...                25-27
        'K_{RM,CI}^a', 'K_{RM,CI}^r', 'K_{RM,Cro}', 'K_{RE}', ...          28-31
        'K_{cro,Cro}', 'K_{cro,CI}', 'K_{cII,Cro}', 'K_{cII,CI}', ...      32-35
        'K_{\lambda,Cro}', 'K_{\lambda,CI}', 'K_{CII}', ...                36-38
        '\deltat', '\tau_{\lambda}'};                                     %39-40

row2c = 1:44;       %Only look at the best 44 fits

%Cook the data to perform all the normalizations
data_ = zeros(max(row2c), 40);
desc_ = cell(40, 1);

for ipar = 1:40
    if ismember(ipar, [4 28 29 33 35 37])    %normalized by Kcro_CI (33)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 33);
        desc_{ipar} = [text{ipar} '/' text{33}];
    elseif ismember(ipar, [6 30 32 34 36])   %normalized by KPRM_Cro (30)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 30);
        desc_{ipar} = [text{ipar} '/' text{30}];
    elseif ismember(ipar, [8 31 38])         %normalized by KPRE (31)
        data_(:, ipar) = data(row2c, ipar) ./ data(row2c, 31);
        desc_{ipar} = [text{ipar} '/' text{31}];
    else                                     %no normalization
        data_(:, ipar) = data(row2c, ipar);
        desc_{ipar} = text{ipar};
    end
end

% Read the decision thresholds
Lysogenic_thresh = zeros(numel(row2c), 1);
Lytic_thresh = zeros(numel(row2c), 1);
% Lytic_thresh = 21.593;
% Lysogenic_thresh = 3.1982;
lytic_perc = 0.5;        %0: lower bound; 1: higher bound
lysog_perc = 0.5;

% lytic_perc = (21.593 - sol{1}.lytic_thresh_min) / (sol{1}.lytic_thresh_max - sol{1}.lytic_thresh_min);
% lysog_perc = (3.1982 - sol{1}.lysog_thresh_min) / (sol{1}.lysog_thresh_max - sol{1}.lysog_thresh_min);

for ii = 1:numel(row2c)
    lytic_thresh = sol{ii}.lytic_thresh_min * (1-lytic_perc) + sol{ii}.lytic_thresh_max * lytic_perc;
    lysog_thresh = sol{ii}.lysog_thresh_min * (1-lysog_perc) + sol{ii}.lysog_thresh_max * lysog_perc;
    
    Lysogenic_thresh(ii) = lysog_thresh;
    Lytic_thresh(ii) = lytic_thresh;
end

%Table S10
lst_tmp = [28 29 30 33 32 35 34 37 36];
desc2_ = cell(numel(lst_tmp), 1);
output2_ = zeros(numel(lst_tmp), 3);

for jparam = 1:numel(lst_tmp)
    
    row_num = lst_tmp(jparam);
    if ismember(row_num, [28 29 33 35 37])
        norm_ = Lysogenic_thresh;
        norm_txt = 'Lysogenic thresh';
    elseif ismember(row_num, [30 32 34 36])
        norm_ = Lytic_thresh;
        norm_txt = 'Lytic thresh';
    end
%     desc_{ipar} = [text{ipar} '/' text{31}];
    desc2_{jparam} = [text{row_num} ' / ' norm_txt];
    output2_(jparam, 1) = round(data_(1, row_num) ./ norm_(1), 2, 'significant');
    
    output2_(jparam, 2) = round(mean(data_(:, row_num) ./ norm_), 2, 'significant');
    output2_(jparam, 3) = round(std(data_(:, row_num) ./ norm_), 2, 'significant');
end

%% Seth run simulations
clear sol
sol(1) = struct();

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
ss = 0;

%ICs
MOI = 1:5;
y0_m1 = [0 0 0 0 0 0 1*convFac/V0];
y0_m2 = [0 0 0 0 0 0 2*convFac/V0];
y0_m3 = [0 0 0 0 0 0 3*convFac/V0];
y0_m4 = [0 0 0 0 0 0 4*convFac/V0];
y0_m5 = [0 0 0 0 0 0 5*convFac/V0];
y0_m10 = [0 0 0 0 0 0 10*convFac/V0];

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, ...
    'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

for i = 1:length(sol)
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod = data(i, 1:9);
    prod(3) = prod(1)*prod(3);
    degr = zeros(8, 1);
    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    [degr(1), degr(2), degr(3), degr(4), degr(5), degr(6), degr(7), degr(8)] = ...
        deal(data(i, 10), data(i, 11), data(i, 12), data(i, 13), data(i, 14), ...
        data(i, 15), data(i, 16), 0); 
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nPR_Cro, nPR_CI, nDeg_CII
    n = [
        data(i, 17);      %nPRM,CI+
        data(i, 18);      %nPRM,CI-
        data(i, 19);      %nPRM,Cro
        data(i, 20);      %nPRE
        data(i, 21);      %nCro,Cro
        data(i, 22);      %nCro,CI
        data(i, 23);      %nCII,Cro
        data(i, 24);      %nCII,CI
        data(i, 25);      %nM,Cro
        data(i, 26);      %nM,CI
        data(i, 27);      %nDeg,CII
         ];
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = zeros(11, 1);
    [K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11)] = ...
        deal(data(i, 28), data(i, 29), data(i, 30), data(i, 31), data(i, 32), ...
        data(i, 33), data(i, 34), data(i, 35), data(i, 36), data(i, 37), ...
        data(i, 38));
    dt = data(i, end-2); %offset
    tau = data(i, end-1);
    kdil = degr(1);
    sol(i).parameters = data(i, :);
    sol(i).prod = prod;
    sol(i).degr = degr;
    sol(i).n = n;
    sol(i).K = K;
    sol(i).dt = dt;
    sol(i).tau = tau;
    
    V = V0.*exp(tspan.*kdil)';
    
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP, ym1_OP] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1 = [tm1_OP, ym1_OP];
    sol(i).m1Num = [tm1_OP, ym1_OP.*V./convFac];
    %MOI=2
    [tm2_OP, ym2_OP] = ode15s(@fv19, tspan, y0_m2(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2 = [tm2_OP, ym2_OP];
    sol(i).m2Num = [tm2_OP, ym2_OP.*V./convFac];
    %MOI=3
    [tm3_OP, ym3_OP] = ode15s(@fv19, tspan, y0_m3(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3 = [tm3_OP, ym3_OP];
    sol(i).m3Num = [tm3_OP, ym3_OP.*V./convFac];
    %MOI=4
    [tm4_OP, ym4_OP] = ode15s(@fv19, tspan, y0_m4(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4 = [tm4_OP, ym4_OP];
    sol(i).m4Num = [tm4_OP, ym4_OP.*V./convFac];
    %MOI=5
    [tm5_OP, ym5_OP] = ode15s(@fv19, tspan, y0_m5(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5 = [tm5_OP, ym5_OP];
    sol(i).m5Num = [tm5_OP, ym5_OP.*V./convFac];
    
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt, ym1_wt] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m1_wt = [tm1_wt, ym1_wt];
    sol(i).m1Num_wt = [tm1_wt, ym1_wt.*V./convFac];
    %MOI=2
    [tm2_wt, ym2_wt] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m2_wt = [tm2_wt, ym2_wt];
    sol(i).m2Num_wt = [tm2_wt, ym2_wt.*V./convFac];
    %MOI=3
    [tm3_wt, ym3_wt] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m3_wt = [tm3_wt, ym3_wt];
    sol(i).m3Num_wt = [tm3_wt, ym3_wt.*V./convFac];
    %MOI=4
    [tm4_wt, ym4_wt] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m4_wt = [tm4_wt, ym4_wt];
    sol(i).m4Num_wt = [tm4_wt, ym4_wt.*V./convFac];
    %MOI=5
    [tm5_wt, ym5_wt] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m5_wt = [tm5_wt, ym5_wt];
    sol(i).m5Num_wt = [tm5_wt, ym5_wt.*V./convFac];  
    
end

%% TY run simulation (rewrite from Seth's codes)
% row2c = 1:44;
row2c = 1;
sol = struct();

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
ss = 0;

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, 'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

for ii = 1:numel(row2c)
%     sol{ii} = struct();
    
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod = data(ii, 1:9);
    prod(3) = prod(1)*prod(3);

    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    degr = [data(ii, 10:16) 0];
    
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, nM_Cro, nM_CI, nDeg_CII
    n = data(ii, 17:27);
    
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = data(ii, 28:38);
    
    dt = data(ii, end-2); %offset
    tau = data(ii, end-1);
    kdil = degr(1);
    sol.parameters = data(ii, :);
    sol.prod = prod;
    sol.degr = degr;
    sol.n = n;
    sol.K = K;
    sol.dt = dt;
    sol.tau = tau;
    
    V = V0.*exp(tspan.*kdil)';
    
    for imoi = 1:5
        y0_ = [zeros(1, 6) imoi * convFac / V0];
        %-------------------- P- ------------------------------------------
        name_1 = ['m' num2str(imoi)];
        name_2 = ['m' num2str(imoi) 'Num'];
        
        [tm_OP, ym_OP] = ode15s(@fv19, tspan, y0_(1:6), options, n([1:8, end]), ...
            prod(1:end-1), degr, K([1:8, end]), imoi, V0, convFac, ss);
        sol.(name_1) = [tm_OP, ym_OP];
        sol.(name_2) = [tm_OP, ym_OP .* V ./ convFac];
        
        %-------------------- P+ ------------------------------------------
        name_3 = ['m' num2str(imoi) '_wt'];
        name_4 = ['m' num2str(imoi) 'Num_wt'];
        
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, tspan, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        sol.(name_3) = [tm_wt, ym_wt];
        sol.(name_4) = [tm_wt, ym_wt .* V ./ convFac];
        
    end
end

%% TY simulation using alternative parameters
clear
sol(1) = struct();

%Normalized data
% CINorm = 56.4071; %KCro_CI
% CroNorm = 52.9038; %KPRM_Cro
% CIINorm = 221.9365; %KPRE
CINorm = 100; %KCro_CI
CroNorm = 10; %KPRM_Cro
CIINorm = 1; %KPRE

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
ss = 0;

%ICs
MOI = 1:5;
y0_m1 = [0 0 0 0 0 0 1*convFac/V0];
y0_m2 = [0 0 0 0 0 0 2*convFac/V0];
y0_m3 = [0 0 0 0 0 0 3*convFac/V0];
y0_m4 = [0 0 0 0 0 0 4*convFac/V0];
y0_m5 = [0 0 0 0 0 0 5*convFac/V0];
y0_m10 = [0 0 0 0 0 0 10*convFac/V0];

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, 'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

for i = 1:length(sol)
    %ASSIGN PARAMS --------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
%     prod = [0.049, 10, 20*0.049, 0.12*0.078*CroNorm, 2.2, 0.078*CroNorm, 5.2, 1.0*0.078*CroNorm, 0.20];
    prod = [0.049, 10, 20*0.049, 0.0092*CINorm, 2.2, 0.078*CroNorm, 5.2, 0.019*CIINorm, 0.20];
%     prod = data(i, 1:9);
%     prod(3) = prod(1)*prod(3);
    
    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    degr = [0.016, 0.10, 0.0083, 0.12, 0.017, 0.10, 1, 0];
%     [degr(1), degr(2), degr(3), degr(4), degr(5), degr(6), degr(7), degr(8)] = ...
%         deal(data(i, 10), data(i, 11), data(i, 12), data(i, 13), data(i, 14), ...
%         data(i, 15), data(i, 16), 0); 
    
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_cro, nCII_CI, nM_Cro, nM_CI, nDeg_CII
    n = [4.0, 6.0, 3.3, 4, 3.1, 4.2, 4, 3.9, 6, 6.0, 1.0];
%     n = [
%         data(i, 17);      %nPRM,CI+
%         data(i, 18);      %nPRM,CI-
%         data(i, 19);      %nPRM,Cro
%         data(i, 20);      %nPRE
%         data(i, 21);      %nCro,Cro
%         data(i, 22);      %nCro,CI
%         data(i, 23);      %nCII,Cro
%         data(i, 24);      %nCII,CI
%         data(i, 25);      %nM,Cro
%         data(i, 26);      %nM,CI
%         data(i, 27);      %nDeg,CII
%          ];
    
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = [1.6*CINorm, 3.1*CINorm, CroNorm, CIINorm, 4.9*CroNorm, CINorm, ...
        4.6*CroNorm, 0.69*CINorm, 19*CroNorm, 3.0*CINorm, 0.84*CIINorm];
%     K = zeros(11, 1);
%     [K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9), K(10), K(11)] = ...
%         deal(data(i, 28), data(i, 29), data(i, 30), data(i, 31), data(i, 32), ...
%         data(i, 33), data(i, 34), data(i, 35), data(i, 36), data(i, 37), ...
%         data(i, 38));
    
    dt = 5;
    tau = 8;
%     dt = data(i, end-2); %offset
%     tau = data(i, end-1);
    kdil = degr(1);
    sol(i).parameters = [prod degr(1:end-1) n K dt tau];
%     sol(i).parameters = data(i, :);
    sol(i).prod = prod;
    sol(i).degr = degr;
    sol(i).n = n;
    sol(i).K = K;
    sol(i).dt = dt;
    sol(i).tau = tau;
    
    V = V0.*exp(tspan.*kdil)';
    
    %P- -------------------------------------------------------------------
    %MOI=1
    [tm1_OP, ym1_OP] = ode15s(@fv19, tspan, y0_m1(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 1, V0, convFac, ss);
    sol(i).m1 = [tm1_OP, ym1_OP];
    sol(i).m1Num = [tm1_OP, ym1_OP.*V./convFac];
    %MOI=2
    [tm2_OP, ym2_OP] = ode15s(@fv19, tspan, y0_m2(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 2, V0, convFac, ss);
    sol(i).m2 = [tm2_OP, ym2_OP];
    sol(i).m2Num = [tm2_OP, ym2_OP.*V./convFac];
    %MOI=3
    [tm3_OP, ym3_OP] = ode15s(@fv19, tspan, y0_m3(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 3, V0, convFac, ss);
    sol(i).m3 = [tm3_OP, ym3_OP];
    sol(i).m3Num = [tm3_OP, ym3_OP.*V./convFac];
    %MOI=4
    [tm4_OP, ym4_OP] = ode15s(@fv19, tspan, y0_m4(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 4, V0, convFac, ss);
    sol(i).m4 = [tm4_OP, ym4_OP];
    sol(i).m4Num = [tm4_OP, ym4_OP.*V./convFac];
    %MOI=5
    [tm5_OP, ym5_OP] = ode15s(@fv19, tspan, y0_m5(1:6), options, n([1:8, end]), ...
        prod(1:end-1), degr, K([1:8, end]), 5, V0, convFac, ss);
    sol(i).m5 = [tm5_OP, ym5_OP];
    sol(i).m5Num = [tm5_OP, ym5_OP.*V./convFac];
    
    
    %WT -------------------------------------------------------------------
    %MOI=1
    [tm1_wt, ym1_wt] = ode15s(@fv19_repv3, tspan, y0_m1, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m1_wt = [tm1_wt, ym1_wt];
    sol(i).m1Num_wt = [tm1_wt, ym1_wt.*V./convFac];
    %MOI=2
    [tm2_wt, ym2_wt] = ode15s(@fv19_repv3, tspan, y0_m2, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m2_wt = [tm2_wt, ym2_wt];
    sol(i).m2Num_wt = [tm2_wt, ym2_wt.*V./convFac];
    %MOI=3
    [tm3_wt, ym3_wt] = ode15s(@fv19_repv3, tspan, y0_m3, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m3_wt = [tm3_wt, ym3_wt];
    sol(i).m3Num_wt = [tm3_wt, ym3_wt.*V./convFac];
    %MOI=4
    [tm4_wt, ym4_wt] = ode15s(@fv19_repv3, tspan, y0_m4, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m4_wt = [tm4_wt, ym4_wt];
    sol(i).m4Num_wt = [tm4_wt, ym4_wt.*V./convFac];
    %MOI=5
    [tm5_wt, ym5_wt] = ode15s(@fv19_repv3, tspan, y0_m5, options, n, prod, ...
        degr, K, tau, V0, convFac);
    sol(i).m5_wt = [tm5_wt, ym5_wt];
    sol(i).m5Num_wt = [tm5_wt, ym5_wt.*V./convFac];  
    
end

%% Seth plot Fig3
%CI, Cro, CII normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

cm = [
    168 222 233;
    101 185 222;
    75 123 190;
    88 78 160;
    7 24 50;
]./255;

fig3_new = figure('Name', 'Phase plane traj.');
dim = [.2 .5 .3 .3];
str_lett = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
numLet = 8;
for i = 1:numLet
    a_lett = annotation('textbox',dim,'String',str_lett{i},'FitBoxToText','on', ...
        'Color', 'k', 'EdgeColor', 'none', ...
        'BackgroundColor', 'none', 'Fontsize', 16);
end
%Phase plane traj ---------------------------------------------------------
subplot(4, 2, 7); %P-   
    %Plot best solution
    p1 = plot(sol(1).m1(:, 5)./CINorm, sol(1).m1(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2(:, 5)./CINorm, sol(1).m2(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3(:, 5)./CINorm, sol(1).m3(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4(:, 5)./CINorm, sol(1).m4(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5(:, 5)./CINorm, sol(1).m5(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    %set(legend, 'NumColumns', 2);
    axis square
subplot(4, 2, 8); %WT
    %Plot best solution
    p1 = plot(sol(1).m1_wt(:, 5)./CINorm, sol(1).m1_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(1, :), 'LineWidth', 1.5); hold on;
    p2 = plot(sol(1).m2_wt(:, 5)./CINorm, sol(1).m2_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(2, :), 'LineWidth', 1.5);
    p3 = plot(sol(1).m3_wt(:, 5)./CINorm, sol(1).m3_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(3, :), 'LineWidth', 1.5);
    p4 = plot(sol(1).m4_wt(:, 5)./CINorm, sol(1).m4_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(4, :), 'LineWidth', 1.5);
    p5 = plot(sol(1).m5_wt(:, 5)./CINorm, sol(1).m5_wt(:, 6)./CroNorm, ...
        '-', 'Color', cm(5, :), 'LineWidth', 1.5);
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
%     set(legend, 'NumColumns', 2);
    xlabel('CI concentration (normalized)');
    ylabel('Cro concentration (normalized)');
    axis square
%CI traj-------------------------------------------------------------------
subplot(4, 2, 1); %P-
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 5)/CINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2(:, 1), sol(1).m2(:, 5)/CINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3(:, 1), sol(1).m3(:, 5)/CINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4(:, 1), sol(1).m4(:, 5)/CINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5(:, 1), sol(1).m5(:, 5)/CINorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('CI concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P-)', 'MOI = 2 (P-)', 'MOI = 3 (P-)', ...
        'MOI = 4 (P-)', 'MOI = 5 (P-)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(4, 2, 2); %P+
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 5)/CINorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 5)/CINorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).m3_wt(:, 5)/CINorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).m4_wt(:, 5)/CINorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).m5_wt(:, 5)/CINorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('CI concentration (normalized)');
    legend([p1, p2, p3, p4, p5], 'MOI = 1 (P+)', 'MOI = 2 (P+)', 'MOI = 3 (P+)', ...
        'MOI = 4 (P+)', 'MOI = 5 (P+)', 'Color', 'none', 'EdgeColor', 'none');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
%Cro traj.-----------------------------------------------------------------
subplot(4, 2, 3); %P-
    p1 = plot(sol(1).m1(:, 1), sol(1).m1(:, 6)/CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2(:, 1), sol(1).m2(:, 6)/CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3(:, 1), sol(1).m3(:, 6)/CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4(:, 1), sol(1).m4(:, 6)/CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5(:, 1), sol(1).m5(:, 6)/CroNorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Cro concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(4, 2, 4); %P+
    p1 = plot(sol(1).m1_wt(:, 1), sol(1).m1_wt(:, 6)/CroNorm, '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2_wt(:, 1), sol(1).m2_wt(:, 6)/CroNorm, '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3_wt(:, 1), sol(1).m3_wt(:, 6)/CroNorm, '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4_wt(:, 1), sol(1).m4_wt(:, 6)/CroNorm, '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5_wt(:, 1), sol(1).m5_wt(:, 6)/CroNorm, '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Cro concentration (normalized)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
%Viral copy number---------------------------------------------------------
subplot(4, 2, 5); %P-
    p1 = plot([0, 65], [1, 1], '-', 'Color', cm(1, :)); hold on;
    p1 = plot([0, 65], [2, 2], '-', 'Color', cm(2, :));
    p1 = plot([0, 65], [3, 3], '-', 'Color', cm(3, :));
    p1 = plot([0, 65], [4, 4], '-', 'Color', cm(4, :));
    p1 = plot([0, 65], [5, 5], '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);
subplot(4, 2, 6); %P+
    p1 = plot(sol(1).m1Num_wt(:, 1), sol(1).m1Num_wt(:, end), '-', 'Color', cm(1, :)); hold on;
    p2 = plot(sol(1).m2Num_wt(:, 1), sol(1).m2Num_wt(:, end), '-', 'Color', cm(2, :));
    p3 = plot(sol(1).m3Num_wt(:, 1), sol(1).m3Num_wt(:, end), '-', 'Color', cm(3, :));
    p4 = plot(sol(1).m4Num_wt(:, 1), sol(1).m4Num_wt(:, end), '-', 'Color', cm(4, :));
    p5 = plot(sol(1).m5Num_wt(:, 1), sol(1).m5Num_wt(:, end), '-', 'Color', cm(5, :));
    xlabel('Time (min)');
    ylabel('Viral copy (#)');
    xlim([0, 65]);
    pbaspect([2, 1, 1]);

%% TY replot Fig3 (CI, Cro, Viral_Copy, Cro_vs_CI) - success
%CI, Cro, CII normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;

fig3_new = figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.45 0.89], ...
    'name', 'Fig 3: Phase plane traj.');

for imoi = 1:5
    %CI trajectory --------------------------------------------------------
%     ylim_CI = 13;
    ylim_CI = 14;
    subplot(4, 2, 1); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 2); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    %Cro trajectory -------------------------------------------------------
    ylim_Cro = 29;
    subplot(4, 2, 3); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 4); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    %Viral copy number ----------------------------------------------------
    ylim_lambda = 350;
    subplot(4, 2, 5); hold on;      %P-
    p_ = plot([0 65], imoi*[1 1], '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 6); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) 'Num_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, end), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
    %Phase plane trajectory -----------------------------------------------
    subplot(4, 2, 7); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 5)/CINorm, data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('CI concentration (normalized)');
        ylabel('Cro concentration (normalized)');
        xlim([0 ylim_CI]);
        ylim([0 ylim_Cro]);
%         lg = legend;
%         set(lg, 'location', 'northwest');
        axis square
    end
    
    subplot(4, 2, 8); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 5)/CINorm, data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('CI concentration (normalized)');
        ylabel('Cro concentration (normalized)');
        xlim([0 ylim_CI]);
        ylim([0 ylim_Cro]);
%         lg = legend;
%         set(lg, 'location', 'northwest');
        axis square
    end
end



%% TY plot new Fig3A-3B (CI, Cro, CII, Viral_copy)
%CI, Cro, CII normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;

fig3_new = figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.45 0.89], ...
    'name', 'Fig 3AB: Phase plane traj.');
ylim_CI = 13;
ylim_Cro = 29;
ylim_CII = 7;
ylim_lambda = 350;
lysog_thresh = 3.198;
lytic_thresh = 21.59;

for imoi = 1:5
    %CI trajectory --------------------------------------------------------
    subplot(4, 2, 1); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        h_ = plot([0 65], lysog_thresh*[1 1], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 2); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        h_ = plot([0 65], lysog_thresh*[1 1], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    %Cro trajectory -------------------------------------------------------
    subplot(4, 2, 3); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        h_ = plot([0 65], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 4); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        h_ = plot([0 65], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    %CII trajectory -------------------------------------------------------
    subplot(4, 2, 5); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi)]);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 7)/CIINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        h_ = plot([0 65], [1 1], '--', 'color', 'm', 'linewidth', 1.5);
        h2_ = plot([0 65], [0.58 0.58], ':', 'color', 'm', 'linewidth', 1.5);
        xlabel('Time (min)');
        ylabel('CII concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CII]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 6); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 7)/CIINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        h_ = plot([0 65], [1 1], '--', 'color', 'm', 'linewidth', 1.5);
        h2_ = plot([0 65], [0.58 0.58], ':', 'color', 'm', 'linewidth', 1.5);
        xlabel('Time (min)');
        ylabel('CII concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CII]);
        pbaspect([2, 1, 1]);
    end
    
    %Viral copy number ----------------------------------------------------
    subplot(4, 2, 7); hold on;      %P-
    p_ = plot([0 65], imoi*[1 1], '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
    subplot(4, 2, 8); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) 'Num_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, end), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
end

%% TY quality control (cI, cro, cII) - Fig2C & FigS6? 
cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;

figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.45 0.7], ...
    'name', 'Supp: RNA traj.');

for imoi = 1:5
    %cI trajectory --------------------------------------------------------
    ylim_cI = 80;
    subplot(3, 2, 1); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi) 'Num']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 2), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cI mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cI]);
        lg = legend;
        set(lg, 'location', 'northeast');
        pbaspect([2, 1, 1]);
    end
    
    subplot(3, 2, 2); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) 'Num_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 2), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cI mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cI]);
        lg = legend;
        set(lg, 'location', 'northeast');
        pbaspect([2, 1, 1]);
    end
    
    %Cro trajectory -------------------------------------------------------
    ylim_cro = 29;
    subplot(3, 2, 3); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi) 'Num']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 3), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cro mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cro]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(3, 2, 4); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) 'Num_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 3), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cro mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cro]);
        pbaspect([2, 1, 1]);
    end
    
    %CII trajectory -------------------------------------------------------
    ylim_cII = 60;
    subplot(3, 2, 5); hold on;      %P-
    data_tmp = sol(1).(['m' num2str(imoi) 'Num']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 4), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P-)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cII mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cII]);
        pbaspect([2, 1, 1]);
    end
    
    subplot(3, 2, 6); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) 'Num_wt']);
    p_ = plot(data_tmp(:, 1), data_tmp(:, 4), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('cII mRNA per cell');
        xlim([0 65]);
        ylim([0 ylim_cII]);
        pbaspect([2, 1, 1]);
    end
    
end

%% TY replot Fig4A
%CI, Cro, CII normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;

fig4_new = figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.45 0.89], ...
    'name', 'Fig 4: PRE activity vs MOI');
yplot = zeros(6, 651);

for imoi = 1:5
    %PRE activity (P+) ----------------------------------------------------
    ylim_PRE = 1;
    subplot(4, 2, 1); hold on;      %P+
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    xplot = data_tmp(:, 1);
    yplot(imoi, :) = 1 ./ (1 + (data_tmp(:, 7)/CIINorm) .^ -sol(1).n(4));   %f_cI,RE([CII])
    
    p_ = plot(xplot, yplot(imoi, :), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        yplot(6, :) = mean(yplot(2:5, :));
        %optional
        p_ = plot(xplot, yplot(6, :), '-', 'color', 'm', ...
            'linewidth', 1.5, 'displayname', 'Average MOI 2-5 (P+)');
        
        xlabel('Time (min)');
        ylabel('P_{RE} activity');
        xlim([0 65]);
        ylim([0 ylim_PRE]);
        pbaspect([2, 1, 1]);
    end
end


%% TY for interest: simulation using alternative PRE activity assumption
%Only simulate best fit, only simulate P+

% row2c = 1:44;
row2c = 1;
sol = struct();

%Constants
convFac = (1e9*1e15)/(6.022e23); %# of um^3 -> nM
V0 = 1; %um^3
ss = 0;

%Solve ODEs
options = odeset('Nonnegative', [], 'RelTol', 1e-6, 'AbsTol', 1e-6);

%tspan
tspan = 0:0.1:65;

%PRE activity
load('Tmp_data.mat', 'PRE_average');

for ii = 1:numel(row2c)
%     sol{ii} = struct();
    
    %ASSIGN PARAMS---------------------------------------------------------
    %cI PRM basal prod, ratio of cI PRM active / basal, ratio cI PRE/PRM basal prod
    %CI translation rate, cro prod rate, Cro translation rate, cII prod rate
    %CII translation rate, replication rate
    prod = data(ii, 1:9);
    prod(3) = prod(1)*prod(3);

    %kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM
    degr = [data(ii, 10:16) 0];
    
    %nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, nM_Cro, nM_CI, nDeg_CII
    n = data(ii, 17:27);
    
    %KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, 
    %KCII_Cro, KCII_CI, KM_Cro, KM_CI, KDeg_CII
    K = data(ii, 28:38);
    
    dt = data(ii, end-2); %offset
    tau = data(ii, end-1);
    kdil = degr(1);
    sol.parameters = data(ii, :);
    sol.prod = prod;
    sol.degr = degr;
    sol.n = n;
    sol.K = K;
    sol.dt = dt;
    sol.tau = tau;
    
    V = V0.*exp(tspan.*kdil)';
    
    for imoi = 1:5
        y0_ = [zeros(1, 6) imoi * convFac / V0];
        
        %-------------------- P+ original ---------------------------------
        name_1 = ['m' num2str(imoi) '_wt'];
        name_2 = ['m' num2str(imoi) 'Num_wt'];
        
        [tm_wt, ym_wt] = ode15s(@fv19_repv3, tspan, y0_, options, n, prod, ...
            degr, K, tau, V0, convFac);
        sol.(name_1) = [tm_wt, ym_wt];
        sol.(name_2) = [tm_wt, ym_wt .* V ./ convFac];
        
        %-------------------- P+ squared window ---------------------------
        name_3 = ['m' num2str(imoi) '_wt_sqr'];
        name_4 = ['m' num2str(imoi) 'Num_wt_sqr'];
        
        %       MOIs:    1    2    3    4    5   Avg
        PRE_window = [ 4.2  2.6  2.0  1.7  1.5  2.4; ...    turn on
                      23.0 22.9 24.8 26.7 28.5 25.18; ...   turn off
                      0.82    1    1    1    1    1];%     amplitude
        t1t2A = PRE_window(:, imoi);
        [tm_wt, ym_wt] = ode15s(@fv19_repsp1, tspan, y0_, options, n, prod, ...
            degr, K, tau, V0, t1t2A, convFac);
        sol.(name_3) = [tm_wt, ym_wt];
        sol.(name_4) = [tm_wt, ym_wt .* V ./ convFac];
        
        %-------------------- P+ scaled (only for MOI = 1) ----------------
        name_5 = ['m' num2str(imoi) '_wt_scl'];
        name_6 = ['m' num2str(imoi) 'Num_wt_scl'];
        
        if imoi == 1
            [tm_wt, ym_wt] = ode15s(@fv19_repsp2, tspan, y0_, options, n, prod, ...
                degr, K, tau, V0, convFac);
            sol.(name_5) = [tm_wt, ym_wt];
            sol.(name_6) = [tm_wt, ym_wt .* V ./ convFac];
        end
        
    end
end


%% TY for interest: plot Fig4A and Fig3
%CI, Cro, CII normalizations
CINorm = sol(1).K(6); %KCro_CI
CroNorm = sol(1).K(3); %KPRM_Cro
CIINorm = sol(1).K(4); %KPRE

cm = [168 222 233;
      101 185 222;
       75 123 190;
       88  78 160;
        7  24  50]./255;

fig4_new = figure('numbertitle', 'off', 'unit', 'normalized', 'position', [0 0 0.98 0.92], ...
    'name', 'Artificial PRE activity curves');
ylim_PRE = 1;
ylim_CI = 14;
ylim_Cro = 30;
ylim_lambda = 350;
lysog_thresh = 3.198;
lytic_thresh = 21.59;

for imoi = 1:5
    %-------------------- P+ original -------------------------------------
    data_tmp = sol(1).(['m' num2str(imoi) '_wt']);
    data_tmp_num = sol(1).(['m' num2str(imoi) 'Num_wt']);
    
    %PRE activity curves
    subplot(3, 5, 1); hold on;
    xplot = data_tmp(:, 1);
    yplot = 1 ./ (1 + (data_tmp(:, 7)/CIINorm) .^ -sol(1).n(4));   %f_cI,RE([CII])
    
    p_ = plot(xplot, yplot, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('P_{RE} activity');
        title('Original fit: main Fig4A, 3B & 3D');
        xlim([0 65]);
        ylim([0 ylim_PRE]);
        pbaspect([2, 1, 1]);
    end
    
    %CI trajectory
    subplot(3, 5, 2); hold on;
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t_ = plot([0 65], lysog_thresh*[1 1], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    %Cro trajectory
    subplot(3, 5, 3); hold on;
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t_ = plot([0 65], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    %Viral copy number
    subplot(3, 5, 4); hold on;
    p_ = plot(data_tmp_num(:, 1), data_tmp_num(:, end), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
    %Phase plane trajectory
    subplot(3, 5, 5); hold on;
    p_ = plot(data_tmp(:, 5)/CINorm, data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t1_ = plot(lysog_thresh*[1 1], [0 ylim_Cro], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        t2_ = plot([0 ylim_CI], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('CI concentration (normalized)');
        ylabel('Cro concentration (normalized)');
        xlim([0 ylim_CI]);
        ylim([0 ylim_Cro]);
        axis square
    end
    
    
    
    %-------------------- P+ squared window -------------------------------
    data_tmp = sol(1).(['m' num2str(imoi) '_wt_sqr']);
    data_tmp_num = sol(1).(['m' num2str(imoi) 'Num_wt_sqr']);
    
    %PRE activity curves
    subplot(3, 5, 6); hold on;
    xplot = data_tmp(:, 1);
    %       MOIs:    1    2    3    4    5   Avg
    PRE_window = [ 4.2  2.6  2.0  1.7  1.5  2.4; ...    turn on
                  23.0 22.9 24.8 26.7 28.5 25.18; ...   turn off
                  0.82    1    1    1    1    1];%     amplitude
    t1t2A = PRE_window(:, imoi);
    yplot = activeWindow(xplot, t1t2A(1), t1t2A(2), t1t2A(3));
    
    p_ = plot(xplot, yplot, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('P_{RE} activity');
        title('Artificial P_{RE} activity window');
        xlim([0 65]);
        ylim([0 ylim_PRE]);
        pbaspect([2, 1, 1]);
    end
    
    %CI trajectory
    subplot(3, 5, 7); hold on;
    p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t_ = plot([0 65], lysog_thresh*[1 1], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
    end
    
    %Cro trajectory
    subplot(3, 5, 8); hold on;
    p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t_ = plot([0 65], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
    end
    
    %Viral copy number
    subplot(3, 5, 9); hold on;
    p_ = plot(data_tmp_num(:, 1), data_tmp_num(:, end), '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
    end
    
    %Phase plane trajectory
    subplot(3, 5, 10); hold on;
    p_ = plot(data_tmp(:, 5)/CINorm, data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
        'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
    if imoi == 5
        t1_ = plot(lysog_thresh*[1 1], [0 ylim_Cro], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        t2_ = plot([0 ylim_CI], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('CI concentration (normalized)');
        ylabel('Cro concentration (normalized)');
        xlim([0 ylim_CI]);
        ylim([0 ylim_Cro]);
        axis square
    end
    
    
    
    %-------------------- P+ scaled (only for MOI = 1) --------------------
    if imoi == 1
        data_tmp = sol(1).(['m' num2str(imoi) '_wt_scl']);
        data_tmp_num = sol(1).(['m' num2str(imoi) 'Num_wt_scl']);
        
        %PRE activity curves
        subplot(3, 5, 11); hold on;
        xplot = data_tmp(:, 1);
        yplot = 1 ./ (1 + (data_tmp(:, 7)/CIINorm) .^ -sol(1).n(4));   %f_cI,RE([CII])
        
        p_ = plot(xplot, 1.2*yplot, '-', 'color', cm(imoi, :), ...
            'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        xlabel('Time (min)');
        ylabel('P_{RE} activity');
        title('Scale up the P_{RE} activity of MOI = 1');
        xlim([0 65]);
        ylim([0 ylim_PRE]);
        pbaspect([2, 1, 1]);
        
        %CI trajectory
        subplot(3, 5, 12); hold on;
        p_ = plot(data_tmp(:, 1), data_tmp(:, 5)/CINorm, '-', 'color', cm(imoi, :), ...
            'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        t_ = plot([0 65], lysog_thresh*[1 1], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        xlabel('Time (min)');
        ylabel('CI concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_CI]);
        pbaspect([2, 1, 1]);
        
        %Cro trajectory
        subplot(3, 5, 13); hold on;
        p_ = plot(data_tmp(:, 1), data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
            'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        t_ = plot([0 65], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('Time (min)');
        ylabel('Cro concentration (normalized)');
        xlim([0 65]);
        ylim([0 ylim_Cro]);
        pbaspect([2, 1, 1]);
        
        %Viral copy number
        subplot(3, 5, 14); hold on;
        p_ = plot(data_tmp_num(:, 1), data_tmp_num(:, end), '-', 'color', cm(imoi, :), ...
            'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        xlabel('Time (min)');
        ylabel('Viral copy (#)');
        xlim([0 65]);
        ylim([0.8 ylim_lambda]);
        set(gca, 'yscale', 'log');
        pbaspect([2, 1, 1]);
        
        %Phase plane trajectory
        subplot(3, 5, 15); hold on;
        p_ = plot(data_tmp(:, 5)/CINorm, data_tmp(:, 6)/CroNorm, '-', 'color', cm(imoi, :), ...
            'linewidth', 1.5, 'displayname', ['MOI = ' num2str(imoi) ' (P+)']);
        t1_ = plot(lysog_thresh*[1 1], [0 ylim_Cro], '--', 'color', 'r', ...
            'linewidth', 1.5, 'displayname', 'Lysogenic threshold');
        t2_ = plot([0 ylim_CI], lytic_thresh*[1 1], '--', 'color', 'g', ...
            'linewidth', 1.5, 'displayname', 'Lytic threshold');
        xlabel('CI concentration (normalized)');
        ylabel('Cro concentration (normalized)');
        xlim([0 ylim_CI]);
        ylim([0 ylim_Cro]);
        axis square
    end
end










