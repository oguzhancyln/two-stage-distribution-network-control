%% test02_multirun_gwo_robustness_fullhybrid.m
% Stage 2 test for the IEEE Access paper:
% Multi-run robustness of the GWO-based safety-gated hybrid strategy.
%
% What it does:
%   - Re-optimizes DNR+CAP windows and PV-active DROOP windows for each run.
%   - Uses the same safety-gate rule for the hybrid case:
%       accept droop only if J_post <= J_pre for that time step.
%   - Reports best/mean/std metrics and runtime across independent runs.
%
% Requirements in MATLAB path/current folder:
%   - MATPOWER with case33bw and case69modified/case69.
%   - 33-bus profiles: myloadspvs.mat with myloads and myPVs/PV.
%   - 69-bus profiles: myloads.mat and mypvs.mat with PV.
%
% Output folder:
%   results_test02_multirun_gwo/
%
% Recommended first debug run:
%   cfg.nRuns = 2; cfg.topoIter = 20; cfg.droopIter = 15;
% Final paper run:
%   cfg.nRuns = 5; cfg.topoIter = 100; cfg.droopIter = 50;

clear; clc; close all;
warning('off','all');

%% ---------------- USER SETTINGS ----------------
cfg.outputDir = 'results_test02_multirun_gwo';
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end

% Final-paper settings. If too slow, reduce nRuns/iterations for debugging.
cfg.nRuns = 5;
cfg.topoAgents = 50;
cfg.topoIter   = 100;
cfg.droopAgents = 40;
cfg.droopIter   = 50;

cfg.verboseIter = true;       % prints GWO progress
cfg.baseSeed = 24001;         % deterministic reproducibility
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambda = 1e5;
cfg.pfmin = 0.90;
cfg.maxCapStep = 5;
cfg.capStepMVAr = 60/1000;   % 60 kVAr per step, in MVAr

cfg.topoWindows  = [1 24; 40 60; 73 96];
cfg.droopWindows = [25 39; 61 72];

systems = makeSystemConfigs();

allRunRows = {};
allAggRows = {};

for sidx = 1:numel(systems)
    sys = systems(sidx);
    fprintf('\n================ Test02 multi-run: %s ================\n', sys.name);
    [mpc0, Lprof, PVprof, defaultOpen] = loadSystemAndProfiles(sys, cfg);
    T = size(Lprof,2);
    nb = size(mpc0.bus,1);
    fprintf('Loaded %s: %d buses, %d branches, default open=%s\n', sys.caseName, nb, size(mpc0.branch,1), mat2str(defaultOpen));
    fprintf('Profiles: load range=[%.3f, %.3f], PV max=%.4f MW/site\n', min(Lprof(:)), max(Lprof(:)), max(PVprof(:)));

    % Reference no-control once
    V_no = simulateNoDroopFullDay(mpc0, Lprof, PVprof, sys, cfg, defaultOpen, zeros(1,numel(sys.capBuses)));
    M_no = computeMetrics(V_no, cfg);

    runMetrics = [];
    convergenceRows = {};

    for r = 1:cfg.nRuns
        seed = cfg.baseSeed + 1000*sidx + r;
        rng(seed,'twister');
        fprintf('\n--- %s run %d/%d | seed=%d ---\n', sys.name, r, cfg.nRuns, seed);
        runTimer = tic;

        openByT = repmat(defaultOpen(:)', T, 1);
        capsByT = zeros(T, numel(sys.capBuses));
        topoBestByWindow = nan(size(cfg.topoWindows,1),1);
        droopAccepted = false(T,1);
        droopPreObj = nan(T,1);
        droopPostObj = nan(T,1);

        % Stage 1: optimize DNR+CAP windows
        for w = 1:size(cfg.topoWindows,1)
            tRange = cfg.topoWindows(w,1):cfg.topoWindows(w,2);
            label = sprintf('%s run %d TOPO+CAP t=%d..%d', sys.name, r, tRange(1), tRange(end));
            [bestOpen, bestCaps, bestObj, conv] = gwoTopologyCap(mpc0, Lprof, PVprof, sys, cfg, tRange, label);
            topoBestByWindow(w) = bestObj;
            openByT(tRange,:) = repmat(bestOpen(:)', numel(tRange), 1);
            capsByT(tRange,:) = repmat(bestCaps(:)', numel(tRange), 1);

            % Carry this topology/capacitor state into the following droop window, if any.
            if w <= size(cfg.droopWindows,1)
                dRange = cfg.droopWindows(w,1):cfg.droopWindows(w,2);
                openByT(dRange,:) = repmat(bestOpen(:)', numel(dRange), 1);
                capsByT(dRange,:) = repmat(bestCaps(:)', numel(dRange), 1);
            end

            for it = 1:numel(conv)
                convergenceRows(end+1,:) = {sys.name, r, sprintf('TOPOCAP_W%d',w), it, conv(it)}; %#ok<SAGROW>
            end
        end

        % DNR+CAP-only voltage matrix
        V_dnrCap = nan(T, nb);
        for t = 1:T
            V_dnrCap(t,:) = runOneTimeStepNoDroop(mpc0, Lprof, PVprof, sys, cfg, t, openByT(t,:), capsByT(t,:));
        end
        M_dnr = computeMetrics(V_dnrCap, cfg);

        % Stage 2: optimize droop at PV-active steps and safety-gate it.
        V_hybrid = V_dnrCap;
        for dw = 1:size(cfg.droopWindows,1)
            for t = cfg.droopWindows(dw,1):cfg.droopWindows(dw,2)
                preV = V_dnrCap(t,:);
                preJ = timeStepObjective(preV, cfg);
                label = sprintf('%s run %d DROOP t=%d', sys.name, r, t);
                [bestEta, bestObj, postV, conv] = gwoDroopOneStep(mpc0, Lprof, PVprof, sys, cfg, t, openByT(t,:), capsByT(t,:), label);
                droopPreObj(t) = preJ;
                droopPostObj(t) = bestObj;
                if bestObj <= preJ
                    V_hybrid(t,:) = postV;
                    droopAccepted(t) = true;
                else
                    V_hybrid(t,:) = preV;
                    droopAccepted(t) = false;
                end
                %#ok<NASGU> bestEta retained for future extension if needed
                for it = 1:numel(conv)
                    convergenceRows(end+1,:) = {sys.name, r, sprintf('DROOP_t%d',t), it, conv(it)}; %#ok<SAGROW>
                end
            end
        end
        runTime = toc(runTimer);
        M_hyb = computeMetrics(V_hybrid, cfg);

        fprintf('%s run %d DONE | J=%.6g | viol=%d | minV=%.5f | maxV=%.5f | accepted droop=%d/%d | runtime=%.1fs\n', ...
            sys.name, r, M_hyb.TotalVoltageDeviation, M_hyb.TotalViolationCount, M_hyb.MinV, M_hyb.MaxV, ...
            nnz(droopAccepted), sum(isInAnyWindow(1:T, cfg.droopWindows)), runTime);

        row = struct();
        row.System = sys.name;
        row.Run = r;
        row.Seed = seed;
        row.NoControl_JV = M_no.TotalVoltageDeviation;
        row.NoControl_Viol = M_no.TotalViolationCount;
        row.DNR_CAP_JV = M_dnr.TotalVoltageDeviation;
        row.DNR_CAP_Viol = M_dnr.TotalViolationCount;
        row.Proposed_JV = M_hyb.TotalVoltageDeviation;
        row.Proposed_Viol = M_hyb.TotalViolationCount;
        row.Proposed_MinV = M_hyb.MinV;
        row.Proposed_MaxV = M_hyb.MaxV;
        row.Proposed_MeanV = M_hyb.MeanV;
        row.DroopAcceptedSteps = nnz(droopAccepted);
        row.RuntimeSeconds = runTime;
        runMetrics = [runMetrics; row]; %#ok<AGROW>

        save(fullfile(cfg.outputDir, sprintf('test02_%s_run%02d.mat', sys.name, r)), ...
            'sys','cfg','seed','openByT','capsByT','V_no','V_dnrCap','V_hybrid', ...
            'M_no','M_dnr','M_hyb','topoBestByWindow','droopAccepted','droopPreObj','droopPostObj');
    end

    % Write per-run CSV
    T_runs = struct2table(runMetrics);
    writetable(T_runs, fullfile(cfg.outputDir, sprintf('test02_multirun_%s_runs.csv', sys.name)));

    % Aggregated robustness statistics
    agg = aggregateRuns(T_runs, sys.name);
    writetable(agg, fullfile(cfg.outputDir, sprintf('test02_multirun_%s_aggregate.csv', sys.name)));

    % Store convergence traces
    if ~isempty(convergenceRows)
        T_conv = cell2table(convergenceRows, 'VariableNames', {'System','Run','Stage','Iteration','BestObjective'});
        writetable(T_conv, fullfile(cfg.outputDir, sprintf('test02_multirun_%s_convergence.csv', sys.name)));
    end

    % Append combined tables
    for i = 1:height(T_runs)
        allRunRows(end+1,:) = table2cell(T_runs(i,:)); %#ok<SAGROW>
    end
    for i = 1:height(agg)
        allAggRows(end+1,:) = table2cell(agg(i,:)); %#ok<SAGROW>
    end
end

% Combined files
if ~isempty(allRunRows)
    varNames = T_runs.Properties.VariableNames;
    T_all_runs = cell2table(allRunRows, 'VariableNames', varNames);
    writetable(T_all_runs, fullfile(cfg.outputDir, 'test02_multirun_all_runs.csv'));
end
if ~isempty(allAggRows)
    aggNames = agg.Properties.VariableNames;
    T_all_agg = cell2table(allAggRows, 'VariableNames', aggNames);
    writetable(T_all_agg, fullfile(cfg.outputDir, 'test02_multirun_aggregate.csv'));
    writeLatexRobustnessTable(T_all_agg, fullfile(cfg.outputDir, 'test02_multirun_aggregate_latex_table.tex'));
end

fprintf('\nDone. Test02 outputs saved under: %s\n', cfg.outputDir);

%% ============================ LOCAL FUNCTIONS ============================
function systems = makeSystemConfigs()
    systems(1).name = '33bus';
    systems(1).caseName = 'case33bw';
    systems(1).fallbackCaseName = 'case33bw';
    systems(1).defaultOpen = [33 34 35 36 37];
    systems(1).pvBuses = [21 24 32 17 18];
    systems(1).capBuses = [18 25 33];
    systems(1).pvRatedMW = 0.20;
    systems(1).loadFile = 'myloadspvs.mat';
    systems(1).pvFile = 'myloadspvs.mat';

    systems(2).name = '69bus';
    systems(2).caseName = 'case69modified';
    systems(2).fallbackCaseName = 'case69';
    systems(2).defaultOpen = [69 70 71 72 73];
    systems(2).pvBuses = [15 18 27 35 46 60];
    systems(2).capBuses = [7 23 30 39 55];
    systems(2).pvRatedMW = 1.00;
    systems(2).loadFile = 'myloads.mat';
    systems(2).pvFile = 'mypvs.mat';
end

function [mpc0, Lprof, PVprof, defaultOpen] = loadSystemAndProfiles(sys, cfg)
    try
        mpc0 = loadcase(sys.caseName);
    catch
        warning('Could not load %s. Trying %s.', sys.caseName, sys.fallbackCaseName);
        mpc0 = loadcase(sys.fallbackCaseName);
    end
    nBr = size(mpc0.branch,1);
    defaultOpen = sys.defaultOpen(sys.defaultOpen <= nBr);
    if numel(defaultOpen) ~= nBr - (size(mpc0.bus,1)-1)
        defaultOpen = find(mpc0.branch(:,11)==0)';
        if isempty(defaultOpen)
            defaultOpen = (nBr-(nBr-(size(mpc0.bus,1)-1))+1):nBr;
        end
    end
    mpc0.branch(:,11) = 1;

    nb = size(mpc0.bus,1); nPV = numel(sys.pvBuses);
    [Lprof, PVprof] = loadProfiles(sys, nb, nPV);

    if strcmpi(sys.name,'33bus')
        mx = max(PVprof(:));
        if mx > 0
            PVprof = PVprof ./ mx .* sys.pvRatedMW;
        end
    end

    if size(Lprof,2) ~= 96 || size(PVprof,2) ~= 96
        error('Profiles must have 96 time steps. Found load %s and PV %s.', mat2str(size(Lprof)), mat2str(size(PVprof)));
    end
    if size(Lprof,1) ~= nb
        error('Load profile bus count mismatch for %s. Expected %d rows, found %d.', sys.name, nb, size(Lprof,1));
    end
    if size(PVprof,1) ~= nPV
        error('PV profile count mismatch for %s. Expected %d rows, found %d.', sys.name, nPV, size(PVprof,1));
    end
end

function [Lprof, PVprof] = loadProfiles(sys, nb, nPV)
    Lprof = [];
    PVprof = [];
    if strcmpi(sys.name,'33bus') && exist(sys.loadFile,'file')
        S = load(sys.loadFile);
        Lprof = getFieldIfExists(S, {'myloads','loads','L','Lprof'});
        PVprof = getFieldIfExists(S, {'myPVs','PV','mypvs','G','PVprof'});
    elseif strcmpi(sys.name,'69bus')
        if exist(sys.loadFile,'file')
            S = load(sys.loadFile);
            Lprof = getFieldIfExists(S, {'myloads','loads','L','Lprof'});
        end
        if exist(sys.pvFile,'file')
            P = load(sys.pvFile);
            PVprof = getFieldIfExists(P, {'PV','myPVs','mypvs','G','PVprof'});
        end
    end

    if isempty(Lprof) || isempty(PVprof)
        warning('Could not load profiles for %s. Generating synthetic profiles for debug only.', sys.name);
        [Lprof, PVprof] = syntheticProfiles(nb, nPV, sys.pvRatedMW);
    end
    Lprof = orientProfile(Lprof, nb, 96);
    PVprof = orientProfile(PVprof, nPV, 96);
end

function X = getFieldIfExists(S, names)
    X = [];
    for i = 1:numel(names)
        if isfield(S, names{i})
            X = S.(names{i});
            return;
        end
    end
end

function X = orientProfile(X, nRows, nCols)
    X = double(X);
    if isequal(size(X), [nRows nCols])
        return;
    elseif isequal(size(X), [nCols nRows])
        X = X';
    else
        error('Profile size %s cannot be oriented to %dx%d.', mat2str(size(X)), nRows, nCols);
    end
end

function [L, PV] = syntheticProfiles(nb, nPV, pvRatedMW)
    T = 96; dt = 0.25; h = (0:T-1)*dt;
    rho = 0.7*ones(1,T);
    for t = 1:T
        x = h(t);
        if x < 5
            rho(t) = 0.7;
        elseif x < 8
            rho(t) = 0.7 + (x-5)/3*(0.9-0.7);
        elseif x < 15
            rho(t) = 0.9;
        elseif x < 19
            rho(t) = 0.9 + (x-15)/4*(1.1-0.9);
        elseif x < 22
            rho(t) = 1.2 - (x-19)/3*(1.1-0.7);
        else
            rho(t) = 0.7;
        end
    end
    L = repmat(rho, nb, 1);
    g = max(0, sin(pi*(h-6)/12));
    g(h<6 | h>18) = 0;
    PV = pvRatedMW * repmat(g, nPV, 1);
end

function Vmat = simulateNoDroopFullDay(mpc0, Lprof, PVprof, sys, cfg, openLines, caps)
    T = size(Lprof,2); nb = size(mpc0.bus,1);
    Vmat = nan(T,nb);
    for t = 1:T
        Vmat(t,:) = runOneTimeStepNoDroop(mpc0, Lprof, PVprof, sys, cfg, t, openLines, caps);
    end
end

function V = runOneTimeStepNoDroop(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps)
    mpc = buildMpcAtTime(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps);
    mpopt = mpoption('verbose',0,'out.all',0);
    r = runpf(mpc, mpopt);
    if ~r.success
        V = 0.5*ones(1,size(mpc0.bus,1));
    else
        V = abs(r.bus(:,8))';
    end
end

function [Vpre, Vpost, Jpre, Jpost] = runOneTimeStepWithDroop(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps, eta)
    mpc = buildMpcAtTime(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps);
    mpopt = mpoption('verbose',0,'out.all',0);
    r0 = runpf(mpc, mpopt);
    if ~r0.success
        Vpre = 0.5*ones(1,size(mpc0.bus,1));
        Vpost = Vpre; Jpre = 1e12; Jpost = 1e12; return;
    end
    Vpre = abs(r0.bus(:,8))';
    Jpre = timeStepObjective(Vpre, cfg);

    nPV = numel(sys.pvBuses);
    Vsp = eta(1:nPV);
    Vband = eta(nPV+1:end);
    baseMVA = mpc.baseMVA;
    Ppv = PVprof(:,t);
    Smax = sys.pvRatedMW / cfg.pfmin;
    if numel(Smax)==1, Smax = Smax*ones(nPV,1); end
    qByInv = zeros(nPV,1);
    for k = 1:nPV
        qCapCircle = sqrt(max(Smax(k)^2 - Ppv(k)^2, 0));
        qCapPF = Ppv(k) * tan(acos(cfg.pfmin));
        qmaxMW = min(qCapCircle, qCapPF);
        QmaxPU = qmaxMW / baseMVA;
        vk = Vpre(sys.pvBuses(k));
        qpu = pvControlSetPointLocal(vk, Vsp(k), Vband(k), QmaxPU);
        qByInv(k) = qpu * baseMVA;
    end
    for k = 1:nPV
        mpc.bus(sys.pvBuses(k),4) = mpc.bus(sys.pvBuses(k),4) - qByInv(k);
    end
    r1 = runpf(mpc, mpopt);
    if ~r1.success
        Vpost = Vpre; Jpost = 1e12;
    else
        Vpost = abs(r1.bus(:,8))';
        Jpost = timeStepObjective(Vpost, cfg);
    end
end

function mpc = buildMpcAtTime(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps)
    mpc = mpc0;
    mpc.branch(:,11) = 1;
    openLines = unique(round(openLines(:)'));
    openLines = openLines(openLines>=1 & openLines<=size(mpc.branch,1));
    mpc.branch(openLines,11) = 0;

    P0 = mpc.bus(:,3); Q0 = mpc.bus(:,4);
    mpc.bus(:,3) = P0 .* Lprof(:,t);
    mpc.bus(:,4) = Q0 .* Lprof(:,t);

    for k = 1:numel(sys.pvBuses)
        mpc.bus(sys.pvBuses(k),3) = mpc.bus(sys.pvBuses(k),3) - PVprof(k,t);
    end
    capSteps = round(capSteps(:)');
    for c = 1:numel(sys.capBuses)
        mpc.bus(sys.capBuses(c),4) = mpc.bus(sys.capBuses(c),4) - capSteps(c)*cfg.capStepMVAr;
    end
end

function [bestOpen, bestCaps, bestObj, conv] = gwoTopologyCap(mpc0, Lprof, PVprof, sys, cfg, tRange, label)
    nb = size(mpc0.bus,1); nBr = size(mpc0.branch,1);
    nOpen = nBr - (nb - 1);
    nCap = numel(sys.capBuses);
    dim = nOpen + nCap;
    lb = [ones(1,nOpen), zeros(1,nCap)];
    ub = [nBr*ones(1,nOpen), cfg.maxCapStep*ones(1,nCap)];

    X = zeros(cfg.topoAgents, dim);
    for a = 1:cfg.topoAgents
        X(a,1:nOpen) = randomRadialOpen(mpc0, nOpen);
        X(a,nOpen+1:end) = randi([0 cfg.maxCapStep], 1, nCap);
    end

    alpha = X(1,:); beta = X(1,:); delta = X(1,:);
    alphaScore = inf; betaScore = inf; deltaScore = inf;
    conv = nan(cfg.topoIter,1);

    for it = 1:cfg.topoIter
        feasible = 0;
        for a = 1:cfg.topoAgents
            [openLines, caps] = decodeTopologyCap(X(a,:), mpc0, nOpen, nCap, cfg);
            X(a,:) = [openLines caps];
            feasible = feasible + 1;
            score = fitnessTopologyCap(mpc0, Lprof, PVprof, sys, cfg, tRange, openLines, caps);
            if score < alphaScore
                deltaScore = betaScore; delta = beta;
                betaScore = alphaScore; beta = alpha;
                alphaScore = score; alpha = X(a,:);
            elseif score < betaScore
                deltaScore = betaScore; delta = beta;
                betaScore = score; beta = X(a,:);
            elseif score < deltaScore
                deltaScore = score; delta = X(a,:);
            end
        end
        aGWO = 2 - it*(2/cfg.topoIter);
        X = updateGWOPositions(X, alpha, beta, delta, aGWO, lb, ub);
        conv(it) = alphaScore;
        if cfg.verboseIter && (mod(it,10)==0 || it==1 || it==cfg.topoIter)
            fprintf('[%s] iter %3d/%3d | best=%12.6g | feasible=%d/%d | open=%s caps=%s\n', ...
                label, it, cfg.topoIter, alphaScore, feasible, cfg.topoAgents, mat2str(alpha(1:nOpen)), mat2str(alpha(nOpen+1:end)));
        end
    end
    [bestOpen, bestCaps] = decodeTopologyCap(alpha, mpc0, nOpen, nCap, cfg);
    bestObj = alphaScore;
end

function score = fitnessTopologyCap(mpc0, Lprof, PVprof, sys, cfg, tRange, openLines, caps)
    score = 0;
    for t = tRange
        V = runOneTimeStepNoDroop(mpc0, Lprof, PVprof, sys, cfg, t, openLines, caps);
        score = score + timeStepObjective(V, cfg);
    end
end

function [bestEta, bestObj, bestV, conv] = gwoDroopOneStep(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps, label)
    nPV = numel(sys.pvBuses);
    dim = 2*nPV;
    lb = [0.95*ones(1,nPV), 0.01*ones(1,nPV)];
    ub = [1.05*ones(1,nPV), 0.05*ones(1,nPV)];
    X = rand(cfg.droopAgents, dim).*(ub-lb) + lb;
    alpha = X(1,:); beta = X(1,:); delta = X(1,:);
    alphaScore = inf; betaScore = inf; deltaScore = inf;
    bestV = nan(1,size(mpc0.bus,1));
    conv = nan(cfg.droopIter,1);
    for it = 1:cfg.droopIter
        feasible = 0;
        for a = 1:cfg.droopAgents
            X(a,:) = max(min(X(a,:), ub), lb);
            [~, Vpost, ~, score] = runOneTimeStepWithDroop(mpc0, Lprof, PVprof, sys, cfg, t, openLines, capSteps, X(a,:));
            if isfinite(score), feasible = feasible + 1; end
            if score < alphaScore
                deltaScore = betaScore; delta = beta;
                betaScore = alphaScore; beta = alpha;
                alphaScore = score; alpha = X(a,:); bestV = Vpost;
            elseif score < betaScore
                deltaScore = betaScore; delta = beta;
                betaScore = score; beta = X(a,:);
            elseif score < deltaScore
                deltaScore = score; delta = X(a,:);
            end
        end
        aGWO = 2 - it*(2/cfg.droopIter);
        X = updateGWOPositions(X, alpha, beta, delta, aGWO, lb, ub);
        conv(it) = alphaScore;
        if cfg.verboseIter && (mod(it,10)==0 || it==1 || it==cfg.droopIter)
            fprintf('[%s] iter %3d/%3d | best=%12.6g | feasible=%d/%d\n', ...
                label, it, cfg.droopIter, alphaScore, feasible, cfg.droopAgents);
        end
    end
    bestEta = alpha;
    bestObj = alphaScore;
end

function Xnew = updateGWOPositions(X, alpha, beta, delta, aGWO, lb, ub)
    Xnew = X;
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            r1 = rand; r2 = rand; A1 = 2*aGWO*r1 - aGWO; C1 = 2*r2;
            Dalpha = abs(C1*alpha(j) - X(i,j)); X1 = alpha(j) - A1*Dalpha;
            r1 = rand; r2 = rand; A2 = 2*aGWO*r1 - aGWO; C2 = 2*r2;
            Dbeta = abs(C2*beta(j) - X(i,j)); X2 = beta(j) - A2*Dbeta;
            r1 = rand; r2 = rand; A3 = 2*aGWO*r1 - aGWO; C3 = 2*r2;
            Ddelta = abs(C3*delta(j) - X(i,j)); X3 = delta(j) - A3*Ddelta;
            Xnew(i,j) = (X1 + X2 + X3)/3;
        end
    end
    Xnew = max(min(Xnew, ub), lb);
end

function [openLines, caps] = decodeTopologyCap(x, mpc0, nOpen, nCap, cfg)
    nBr = size(mpc0.branch,1);
    openLines = unique(round(x(1:nOpen)),'stable');
    openLines = openLines(openLines>=1 & openLines<=nBr);
    if numel(openLines) ~= nOpen || ~isRadialAfterOpening(mpc0, openLines)
        openLines = randomRadialOpen(mpc0, nOpen);
    end
    caps = round(x(nOpen+1:nOpen+nCap));
    caps = max(min(caps, cfg.maxCapStep), 0);
end

function openLines = randomRadialOpen(mpc0, nOpen)
    nBr = size(mpc0.branch,1);
    for tries = 1:10000
        openLines = sort(randperm(nBr, nOpen));
        if isRadialAfterOpening(mpc0, openLines)
            return;
        end
    end
    error('Could not sample a radial opening after many tries. Check case topology.');
end

function ok = isRadialAfterOpening(mpc0, openLines)
    nb = size(mpc0.bus,1); nBr = size(mpc0.branch,1);
    status = ones(nBr,1); status(openLines) = 0;
    closed = mpc0.branch(status==1,1:2);
    if size(closed,1) ~= nb-1
        ok = false; return;
    end
    visited = false(nb,1);
    q = 1; visited(1) = true;
    while ~isempty(q)
        u = q(1); q(1) = [];
        idx = closed(:,1)==u | closed(:,2)==u;
        edges = closed(idx,:);
        for e = 1:size(edges,1)
            v = edges(e,1) + edges(e,2) - u;
            if v>=1 && v<=nb && ~visited(v)
                visited(v) = true; q(end+1) = v; %#ok<AGROW>
            end
        end
    end
    ok = all(visited);
end

function Qcmd = pvControlSetPointLocal(Vmeas, Vsp, Vband, Qmax)
    if Vband <= 0, Vband = 0.01; end
    K = Qmax ./ Vband;
    dV = Vmeas - Vsp;
    Qcmd = -K .* dV;
    if dV <= -Vband
        Qcmd = Qmax;
    elseif dV >= Vband
        Qcmd = -Qmax;
    end
    Qcmd = max(min(Qcmd, Qmax), -Qmax);
end

function J = timeStepObjective(V, cfg)
    under = max(cfg.Vmin - V, 0);
    over  = max(V - cfg.Vmax, 0);
    J = sum((V-1).^2) + cfg.lambda * sum(under.^2 + over.^2);
end

function M = computeMetrics(Vmat, cfg)
    under = Vmat < cfg.Vmin;
    over  = Vmat > cfg.Vmax;
    M.MinV = min(Vmat(:));
    M.MaxV = max(Vmat(:));
    M.MeanV = mean(Vmat(:));
    M.TotalVoltageDeviation = sum((Vmat(:)-1).^2);
    M.UnderVoltageCount = nnz(under);
    M.OverVoltageCount = nnz(over);
    M.TotalViolationCount = nnz(under | over);
end

function mask = isInAnyWindow(tVec, windows)
    mask = false(size(tVec));
    for i = 1:size(windows,1)
        mask = mask | (tVec >= windows(i,1) & tVec <= windows(i,2));
    end
end

function agg = aggregateRuns(T_runs, systemName)
    valsJ = T_runs.Proposed_JV;
    valsViol = T_runs.Proposed_Viol;
    valsMinV = T_runs.Proposed_MinV;
    valsMaxV = T_runs.Proposed_MaxV;
    valsRuntime = T_runs.RuntimeSeconds;
    agg = table({systemName}, height(T_runs), min(valsJ), mean(valsJ), std(valsJ), ...
        min(valsViol), mean(valsViol), std(valsViol), ...
        mean(valsMinV), mean(valsMaxV), mean(valsRuntime), std(valsRuntime), ...
        'VariableNames', {'System','Runs','Best_JV','Mean_JV','Std_JV', ...
        'Best_Viol','Mean_Viol','Std_Viol','Mean_MinV','Mean_MaxV','Mean_Runtime_s','Std_Runtime_s'});
end

function writeLatexRobustnessTable(T, outFile)
    fid = fopen(outFile,'w');
    fprintf(fid, '%% Auto-generated by test02_multirun_gwo_robustness_fullhybrid.m\n');
    fprintf(fid, '\\begin{table}[!t]\n\\centering\n');
    fprintf(fid, '\\caption{Multi-run robustness of the safety-gated GWO hybrid strategy.}\n');
    fprintf(fid, '\\label{tab:gwo_multirun_robustness}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrrr}\n\\hline\n');
    fprintf(fid, 'System & Runs & Best $J_V$ & Mean $J_V$ & Std. $J_V$ & Mean viol. & Mean runtime (s) \\\\ \\hline\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %d & %.4f & %.4f & %.4f & %.1f & %.1f \\\\ \n', ...
            T.System{i}, T.Runs(i), T.Best_JV(i), T.Mean_JV(i), T.Std_JV(i), T.Mean_Viol(i), T.Mean_Runtime_s(i));
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end
