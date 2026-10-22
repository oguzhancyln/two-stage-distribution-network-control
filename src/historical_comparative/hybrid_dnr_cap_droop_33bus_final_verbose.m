%% hybrid_dnr_cap_droop_33bus_final_verbose.m
% Final verbose single-file MATLAB implementation for the IEEE 33-bus study.
%
% Purpose
%   1) Builds/loads 96-step load and PV profiles.
%   2) Runs the no-control baseline.
%   3) Optimizes topology reconfiguration + capacitor switching with GWO
%      over selected time windows.
%   4) Optimizes PV smart-inverter Volt/VAr setpoint and band parameters
%      during PV-active windows.
%   5) Exports voltages, selected controls, metrics, and figures for the paper.
%
% Dependencies
%   - MATPOWER must be on the MATLAB path.
%   - case33bw must be available through loadcase('case33bw').
%   - Optional: myloadspvs.mat in the same folder, containing:
%       myloads : 33 x 96 load multipliers
%       myPVs   : 5 x 96 normalized PV profiles OR MW profiles
%
% Notes
%   - This file replaces the fragmented main1_24, main25_39, main40_60,
%     main61_72, main73_96, humhum, sample_k, isRadial, and pvControlSetPoint
%     workflow with one reproducible script.
%   - By default, myPVs is interpreted as a normalized PV profile when its
%     maximum value is <= 1.2, and is multiplied by cfg.PV_RATED_MW.
%   - To reproduce the older prototype behavior where myPVs was subtracted
%     directly as MW, set cfg.PV_PROFILE_MODE = 'MW'.
%
% Author-ready version: May 2026

clear; clc; close all;
warning('off','all');
rng(42);                         % reproducible final paper run; use rng('shuffle') for random trials

%% ---------------------------- USER SETTINGS ----------------------------
cfg.caseName        = 'case33bw';  % MATPOWER case
cfg.outputDir       = 'results_hybrid33';
cfg.profileFile     = 'myloadspvs.mat';

% Study horizon
cfg.T               = 96;          % 24 h at 15-min intervals
cfg.dtHours         = 0.25;
cfg.Vmin            = 0.95;
cfg.Vmax            = 1.05;
cfg.penaltyWeight   = 1e5;

% PV and inverter placement for IEEE 33-bus case
cfg.pvBuses         = [21 24 32 17 18];
cfg.PV_RATED_MW     = 0.20;        % 200 kW per PV site for normalized profile mode
cfg.PV_PROFILE_MODE = 'auto';      % 'auto', 'normalized', or 'MW'
cfg.PFmin           = 0.90;        % minimum power factor for Q capability
cfg.inverterSOverP  = 1/cfg.PFmin; % S_rating = P_rated / PFmin

% Capacitor banks: 0..capMaxSteps steps, each step capStepMVAr
cfg.capBuses        = [18 25 33];
cfg.capStepMVAr     = 60/1000;     % 60 kVAr = 0.060 MVAr
cfg.capMaxSteps     = 5;

% Operating windows used in the original prototype
cfg.reconfigWindows = {1:24, 40:60, 73:96};     % topology + capacitors
cfg.inverterWindows = {25:39, 61:72};           % inverter droop using nearest previous topology/caps

% FINAL GWO settings for paper-quality run.
% Topology/capacitor stage is harder, therefore it uses a larger population
% and more iterations. Droop stage is run separately for each PV-active time
% step, so 40 x 50 is a good final compromise between quality and runtime.
cfg.gwoTopologyAgents = 50;
cfg.gwoTopologyIter   = 100;
cfg.gwoDroopAgents    = 40;
cfg.gwoDroopIter      = 50;

% Verbose progress reporting. Every GWO iteration prints the current best
% objective, decoded control variables, and elapsed time.
cfg.verboseGWO        = true;
cfg.saveCheckpoint    = true;

% Droop parameter bounds
cfg.VspLB          = 0.95;
cfg.VspUB          = 1.05;
cfg.VbandLB        = 0.01;
cfg.VbandUB        = 0.05;

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

%% ----------------------------- LOAD SYSTEM -----------------------------
mpc0 = loadcase(cfg.caseName);
nb   = size(mpc0.bus, 1);
nbr  = size(mpc0.branch, 1);
nOpen = nbr - (nb - 1);

if nb ~= 33
    error('This clean script is configured for the IEEE 33-bus case. Loaded case has %d buses.', nb);
end
if nOpen <= 0
    error('The case does not have enough branches for reconfiguration. nOpen=%d.', nOpen);
end

fprintf('Loaded %s: %d buses, %d branches, %d branches opened for radial operation.\n', ...
    cfg.caseName, nb, nbr, nOpen);

%% --------------------------- LOAD/BUILD PROFILES -----------------------
[loadMult, pvMW] = loadOrCreateProfiles(cfg, nb, numel(cfg.pvBuses));

%% ---------------------------- BASELINE RUN -----------------------------
base.openLines = [];
base.capSteps  = zeros(1, numel(cfg.capBuses));
base.usePVQ    = false;

V_no = nan(cfg.T, nb);
J_no = nan(cfg.T, 1);
for t = 1:cfg.T
    [V_no(t,:), success] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, base.openLines, base.capSteps, [], []);
    if ~success
        warning('Baseline PF did not converge at time step %d.', t);
    end
    J_no(t) = voltageObjective(V_no(t,:), cfg);
end

%% -------------------- STAGE 1: DNR + CAPACITOR GWO ---------------------
controlled.openLinesByWindow = cell(numel(cfg.reconfigWindows), 1);
controlled.capStepsByWindow  = cell(numel(cfg.reconfigWindows), 1);
controlled.topologyCurve     = cell(numel(cfg.reconfigWindows), 1);
controlled.topologyScore     = nan(numel(cfg.reconfigWindows), 1);

V_hybrid = nan(cfg.T, nb);
controlLog = repmat(struct('mode','','openLines',[],'capSteps',[],'Vsp',[],'Vband',[],'score',nan), cfg.T, 1);

for w = 1:numel(cfg.reconfigWindows)
    timeIdx = cfg.reconfigWindows{w};
    fprintf('\nOptimizing topology + capacitors for window %d: t=%d..%d\n', w, timeIdx(1), timeIdx(end));

    [bestOpen, bestCaps, bestScore, curve] = optimizeTopologyCapGWO(mpc0, cfg, loadMult, pvMW, timeIdx, nOpen);

    controlled.openLinesByWindow{w} = bestOpen;
    controlled.capStepsByWindow{w}  = bestCaps;
    controlled.topologyCurve{w}     = curve;
    controlled.topologyScore(w)     = bestScore;

    fprintf('  Best open lines: %s\n', mat2str(bestOpen));
    fprintf('  Best cap steps : %s\n', mat2str(bestCaps));
    fprintf('  Best objective : %.6g\n', bestScore);

    if cfg.saveCheckpoint
        save(fullfile(cfg.outputDir, sprintf('checkpoint_topology_window_%d.mat', w)), ...
            'cfg', 'w', 'timeIdx', 'bestOpen', 'bestCaps', 'bestScore', 'curve');
    end

    for t = timeIdx
        [V_hybrid(t,:), success] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, bestOpen, bestCaps, [], []);
        if ~success
            warning('Hybrid PF did not converge in reconfiguration window %d, time %d.', w, t);
        end
        controlLog(t).mode      = 'DNR+CAP';
        controlLog(t).openLines = bestOpen;
        controlLog(t).capSteps  = bestCaps;
        controlLog(t).score     = voltageObjective(V_hybrid(t,:), cfg);
    end
end

%% -------------------- STAGE 2: SMART-INVERTER DROOP GWO ----------------
% Window mapping: first inverter window uses reconfiguration window 1;
% second inverter window uses reconfiguration window 2.
inverterParentWindow = [1 2];
allDroopCurves = cell(numel(cfg.inverterWindows), 1);

for w = 1:numel(cfg.inverterWindows)
    timeIdx = cfg.inverterWindows{w};
    pWin = inverterParentWindow(w);
    fixedOpen = controlled.openLinesByWindow{pWin};
    fixedCaps = controlled.capStepsByWindow{pWin};

    fprintf('\nOptimizing inverter droop for window %d: t=%d..%d using topology window %d\n', ...
        w, timeIdx(1), timeIdx(end), pWin);

    droopCurveWindow = nan(numel(timeIdx), cfg.gwoDroopIter);
    row = 0;
    for t = timeIdx
        row = row + 1;
        [bestVsp, bestVband, bestScore, curve, Vbest, VpreDroop] = optimizeDroopGWO(mpc0, cfg, loadMult, pvMW, t, fixedOpen, fixedCaps);
        V_hybrid(t,:) = Vbest;
        droopCurveWindow(row,:) = curve;

        controlLog(t).mode      = 'DROOP';
        controlLog(t).openLines = fixedOpen;
        controlLog(t).capSteps  = fixedCaps;
        controlLog(t).Vsp       = bestVsp;
        controlLog(t).Vband     = bestVband;
        controlLog(t).score     = bestScore;

        fprintf('  t=%02d: obj=%.6g, minV pre=%.4f, minV post=%.4f\n', ...
            t, bestScore, min(VpreDroop), min(Vbest));
    end
    allDroopCurves{w} = droopCurveWindow;
    if cfg.saveCheckpoint
        save(fullfile(cfg.outputDir, sprintf('checkpoint_droop_window_%d.mat', w)), ...
            'cfg', 'w', 'timeIdx', 'droopCurveWindow', 'controlLog', 'V_hybrid');
    end
end

%% -------------------- FILL ANY UNASSIGNED STEPS SAFELY -----------------
% In case windows are later edited and leave gaps, use the no-control value
% for unassigned time steps to keep exported arrays complete.
missing = any(isnan(V_hybrid), 2);
if any(missing)
    fprintf('\nWarning: %d time steps were not assigned by the hybrid schedule; filling with baseline.\n', nnz(missing));
    V_hybrid(missing,:) = V_no(missing,:);
    for t = find(missing)'
        controlLog(t).mode = 'BASELINE-FILL';
    end
end

%% ----------------------------- METRICS ---------------------------------
metricsBase   = computeMetrics(V_no, cfg);
metricsHybrid = computeMetrics(V_hybrid, cfg);

summaryTable = table( ...
    {'No control'; 'Hybrid DNR+CAP+DROOP'}, ...
    [metricsBase.minV; metricsHybrid.minV], ...
    [metricsBase.maxV; metricsHybrid.maxV], ...
    [metricsBase.meanV; metricsHybrid.meanV], ...
    [metricsBase.totalDeviation; metricsHybrid.totalDeviation], ...
    [metricsBase.totalPenalty; metricsHybrid.totalPenalty], ...
    [metricsBase.underVoltageCount; metricsHybrid.underVoltageCount], ...
    [metricsBase.overVoltageCount; metricsHybrid.overVoltageCount], ...
    [metricsBase.violationCount; metricsHybrid.violationCount], ...
    'VariableNames', {'Case','MinV','MaxV','MeanV','TotalVoltageDeviation','TotalViolationPenalty','UnderVoltageCount','OverVoltageCount','TotalViolationCount'});

disp(summaryTable);
writetable(summaryTable, fullfile(cfg.outputDir, 'summary_33bus.csv'));
writeLatexSummaryTable(summaryTable, fullfile(cfg.outputDir, 'summary_33bus_latex_table.tex'));

%% ----------------------------- EXPORTS ---------------------------------
save(fullfile(cfg.outputDir, 'hybrid33_results.mat'), ...
    'cfg', 'V_no', 'V_hybrid', 'controlLog', 'controlled', 'summaryTable', 'loadMult', 'pvMW');

exportControlSchedule(controlLog, cfg, fullfile(cfg.outputDir, 'control_schedule_33bus.csv'));
plotAndSaveFigures(cfg, V_no, V_hybrid, controlled, allDroopCurves, pvMW, loadMult);

fprintf('\nDone. Results saved under: %s\n', cfg.outputDir);
fprintf('Use summary_33bus.csv and summary_33bus_latex_table.tex in the paper.\n');

%% =======================================================================
%%                              FUNCTIONS
%% =======================================================================

function [loadMult, pvMW] = loadOrCreateProfiles(cfg, nb, nPV)
    if exist(cfg.profileFile, 'file') == 2
        S = load(cfg.profileFile);
        if isfield(S, 'myloads') && isfield(S, 'myPVs') && size(S.myloads,1) == nb && size(S.myloads,2) >= cfg.T
            loadMult = S.myloads(:, 1:cfg.T);
            rawPV = S.myPVs(1:nPV, 1:cfg.T);
            switch lower(cfg.PV_PROFILE_MODE)
                case 'normalized'
                    pvMW = cfg.PV_RATED_MW * rawPV;
                case 'mw'
                    pvMW = rawPV;
                otherwise
                    if max(rawPV(:)) <= 1.2
                        pvMW = cfg.PV_RATED_MW * rawPV;
                    else
                        pvMW = rawPV;
                    end
            end
            fprintf('Loaded profiles from %s. PV max = %.3f MW/site.\n', cfg.profileFile, max(pvMW(:)));
            return;
        end
    end

    fprintf('Profile file not found or incompatible. Creating synthetic 96-step profiles.\n');
    hours = (0:cfg.T-1) * cfg.dtHours;
    rho = zeros(1, cfg.T);
    for t = 1:cfg.T
        h = hours(t);
        if h < 5
            rho(t) = 0.7;
        elseif h < 8
            rho(t) = 0.7 + (h-5)/3*(0.9-0.7);
        elseif h < 15
            rho(t) = 0.9;
        elseif h < 19
            rho(t) = 0.9 + (h-15)/4*(1.1-0.9);
        elseif h < 22
            rho(t) = 1.2 - (h-19)/3*(1.1-0.7);
        else
            rho(t) = 0.7;
        end
    end
    loadMult = zeros(nb, cfg.T);
    for i = 1:nb
        scale = 0.9 + 0.2*rand();
        noise = 0.02*randn(1, cfg.T);
        loadMult(i,:) = min(1.1, max(0.7, scale*rho + noise));
    end

    g = max(0, sin(pi*(hours-6)/12));
    g(hours < 6 | hours > 18) = 0;
    pvMW = zeros(nPV, cfg.T);
    for k = 1:nPV
        scale = 0.9 + 0.2*rand();
        noise = 0.02*randn(1, cfg.T);
        prof = min(1.1, max(0, scale*g + noise));
        pvMW(k,:) = cfg.PV_RATED_MW * prof;
    end
end

function [V, success] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, openLines, capSteps, Vsp, Vband)
    PD = 3; QD = 4; VM = 8; BR_STATUS = 11;
    mpc = mpc0;
    mpc.branch(:, BR_STATUS) = 1;
    if ~isempty(openLines)
        mpc.branch(openLines, BR_STATUS) = 0;
    end

    % Apply time-varying load.
    mpc.bus(:, PD) = mpc0.bus(:, PD) .* loadMult(:, t);
    mpc.bus(:, QD) = mpc0.bus(:, QD) .* loadMult(:, t);

    % PV active power is represented as negative demand.
    for k = 1:numel(cfg.pvBuses)
        mpc.bus(cfg.pvBuses(k), PD) = mpc.bus(cfg.pvBuses(k), PD) - pvMW(k,t);
    end

    % Capacitor banks reduce reactive demand.
    if ~isempty(capSteps)
        for c = 1:numel(cfg.capBuses)
            mpc.bus(cfg.capBuses(c), QD) = mpc.bus(cfg.capBuses(c), QD) - capSteps(c)*cfg.capStepMVAr;
        end
    end

    mpopt = mpoption('verbose', 0, 'out.all', 0);

    % If Vsp/Vband is provided, first run a pre-droop PF to get local voltage,
    % then apply a one-shot Volt/VAr command and run PF again.
    if ~isempty(Vsp) && ~isempty(Vband)
        r0 = runpf(mpc, mpopt);
        if ~isfield(r0, 'success') || r0.success == 0
            V = nan(1, size(mpc.bus,1)); success = false; return;
        end
        V0 = abs(r0.bus(:, VM));
        Qmax = dynamicQmaxMVAr(cfg, pvMW(:,t));
        for k = 1:numel(cfg.pvBuses)
            qCmd = droopQcmd(V0(cfg.pvBuses(k)), Vsp(k), Vband(k), Qmax(k));
            mpc.bus(cfg.pvBuses(k), QD) = mpc.bus(cfg.pvBuses(k), QD) - qCmd;
        end
    end

    r = runpf(mpc, mpopt);
    success = isfield(r, 'success') && r.success == 1;
    if success
        V = abs(r.bus(:, VM))';
    else
        V = nan(1, size(mpc.bus,1));
    end
end

function qMax = dynamicQmaxMVAr(cfg, pPV)
    % Apparent rating based on installed PV rating.
    S = cfg.PV_RATED_MW * cfg.inverterSOverP * ones(size(pPV));
    qHeadroom = sqrt(max(S.^2 - pPV.^2, 0));
    qPF = pPV .* tan(acos(cfg.PFmin));
    qMax = min(qHeadroom, qPF);
end

function Qcmd = droopQcmd(Vmeas, Vsp, Vband, Qmax)
    if Qmax <= 0 || Vband <= 0
        Qcmd = 0; return;
    end
    K = Qmax ./ Vband;
    dV = Vmeas - Vsp;
    Qcmd = -K .* dV;
    if dV <= -Vband
        Qcmd = Qmax;
    elseif dV >= Vband
        Qcmd = -Qmax;
    end
    Qcmd = min(Qmax, max(-Qmax, Qcmd));
end

function J = voltageObjective(V, cfg)
    if any(isnan(V)) || any(isinf(V))
        J = 1e12; return;
    end
    base = sum((V - 1).^2);
    under = max(cfg.Vmin - V, 0);
    over  = max(V - cfg.Vmax, 0);
    J = base + cfg.penaltyWeight * sum(under.^2 + over.^2);
end

function [bestOpen, bestCaps, bestScore, curve] = optimizeTopologyCapGWO(mpc0, cfg, loadMult, pvMW, timeIdx, nOpen)
    nb = size(mpc0.bus,1);
    nbr = size(mpc0.branch,1);
    nCaps = numel(cfg.capBuses);
    dim = nOpen + nCaps;

    lb = [ones(1,nOpen), zeros(1,nCaps)];
    ub = [nbr*ones(1,nOpen), cfg.capMaxSteps*ones(1,nCaps)];

    nAgents = cfg.gwoTopologyAgents;
    maxIter = cfg.gwoTopologyIter;
    X = zeros(nAgents, dim);
    for i = 1:nAgents
        open = randomRadialOpenSet(mpc0, nOpen);
        caps = randi([0 cfg.capMaxSteps], 1, nCaps);
        X(i,:) = [open caps];
    end

    progressLabel = sprintf('TOPO+CAP window t=%d..%d', timeIdx(1), timeIdx(end));
    decodeFcn = @(x) topologyProgressDecode(x, mpc0, nOpen, cfg);
    [bestScore, bestX, curve] = runGWO(X, lb, ub, maxIter, @(x) topologyFitness(x), cfg, progressLabel, decodeFcn);
    [bestOpen, bestCaps] = decodeTopologyCandidate(bestX, mpc0, nOpen, cfg);

    function J = topologyFitness(x)
        [openLines, capSteps] = decodeTopologyCandidate(x, mpc0, nOpen, cfg);
        if isempty(openLines) || ~isRadialAfterOpening(mpc0, openLines)
            J = 1e12; return;
        end
        J = 0;
        for tt = timeIdx
            [V, success] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, tt, openLines, capSteps, [], []);
            if ~success
                J = J + 1e10;
            else
                J = J + voltageObjective(V, cfg);
            end
        end
    end
end

function [bestVsp, bestVband, bestScore, curve, Vbest, VpreDroop] = optimizeDroopGWO(mpc0, cfg, loadMult, pvMW, t, openLines, capSteps)
    nPV = numel(cfg.pvBuses);
    dim = 2*nPV;
    lb = [cfg.VspLB*ones(1,nPV), cfg.VbandLB*ones(1,nPV)];
    ub = [cfg.VspUB*ones(1,nPV), cfg.VbandUB*ones(1,nPV)];

    nAgents = cfg.gwoDroopAgents;
    X = zeros(nAgents, dim);
    for i = 1:nAgents
        X(i,:) = lb + rand(1,dim).*(ub-lb);
    end

    progressLabel = sprintf('DROOP time step t=%d', t);
    decodeFcn = @(x) droopProgressDecode(x, nPV);
    [bestScore, bestX, curve] = runGWO(X, lb, ub, cfg.gwoDroopIter, @(x) droopFitness(x), cfg, progressLabel, decodeFcn);
    bestVsp = bestX(1:nPV);
    bestVband = bestX(nPV+1:end);
    [Vbest, ~] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, openLines, capSteps, bestVsp, bestVband);
    [VpreDroop, ~] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, openLines, capSteps, [], []);

    function J = droopFitness(x)
        Vsp = x(1:nPV);
        Vband = x(nPV+1:end);
        [V, success] = runPowerFlowAtTime(mpc0, cfg, loadMult, pvMW, t, openLines, capSteps, Vsp, Vband);
        if ~success
            J = 1e12;
        else
            J = voltageObjective(V, cfg);
        end
    end
end

function [alphaScore, alphaX, curve] = runGWO(X, lb, ub, maxIter, fitnessFcn, cfg, progressLabel, decodeFcn)
    % Grey Wolf Optimizer with per-iteration progress printing.
    %
    % progressLabel : text printed before each iteration line
    % decodeFcn     : function handle returning a compact text summary of alphaX

    if nargin < 6 || isempty(cfg)
        cfg.verboseGWO = true;
    end
    if nargin < 7 || isempty(progressLabel)
        progressLabel = 'GWO';
    end
    if nargin < 8 || isempty(decodeFcn)
        decodeFcn = @(x) '';
    end

    nAgents = size(X,1);
    dim = size(X,2);
    alphaX = zeros(1,dim); alphaScore = inf;
    betaX  = zeros(1,dim); betaScore  = inf;
    deltaX = zeros(1,dim); deltaScore = inf;
    curve = nan(1, maxIter);

    tStart = tic;
    lastBest = inf;

    for iter = 1:maxIter
        nFeasibleLike = 0;
        for i = 1:nAgents
            X(i,:) = reflectBounds(X(i,:), lb, ub);
            fit = fitnessFcn(X(i,:));
            if isfinite(fit) && fit < 1e11
                nFeasibleLike = nFeasibleLike + 1;
            end

            if fit < alphaScore
                deltaScore = betaScore;  deltaX = betaX;
                betaScore  = alphaScore; betaX  = alphaX;
                alphaScore = fit;        alphaX = X(i,:);
            elseif fit < betaScore
                deltaScore = betaScore; deltaX = betaX;
                betaScore = fit;        betaX = X(i,:);
            elseif fit < deltaScore
                deltaScore = fit;       deltaX = X(i,:);
            end
        end

        a = 2 - iter*(2/maxIter);
        for i = 1:nAgents
            for j = 1:dim
                r1 = rand(); r2 = rand();
                A1 = 2*a*r1 - a; C1 = 2*r2;
                Dalpha = abs(C1*alphaX(j) - X(i,j));
                X1 = alphaX(j) - A1*Dalpha;

                r1 = rand(); r2 = rand();
                A2 = 2*a*r1 - a; C2 = 2*r2;
                Dbeta = abs(C2*betaX(j) - X(i,j));
                X2 = betaX(j) - A2*Dbeta;

                r1 = rand(); r2 = rand();
                A3 = 2*a*r1 - a; C3 = 2*r2;
                Ddelta = abs(C3*deltaX(j) - X(i,j));
                X3 = deltaX(j) - A3*Ddelta;

                X(i,j) = (X1 + X2 + X3)/3;
            end
            X(i,:) = reflectBounds(X(i,:), lb, ub);
        end

        curve(iter) = alphaScore;

        if cfg.verboseGWO
            if isinf(lastBest)
                improvementPct = NaN;
            else
                improvementPct = 100*(lastBest - alphaScore)/max(abs(lastBest), eps);
            end
            lastBest = alphaScore;
            alphaTxt = decodeFcn(alphaX);
            fprintf('[%s] iter %3d/%3d | best=%12.6g | beta=%12.6g | delta=%12.6g | feasible=%2d/%2d | improv=%8.3f%% | elapsed=%7.1fs | %s\n', ...
                progressLabel, iter, maxIter, alphaScore, betaScore, deltaScore, nFeasibleLike, nAgents, improvementPct, toc(tStart), alphaTxt);
        end
    end
end

function txt = topologyProgressDecode(x, mpc, nOpen, cfg)
    try
        [openLines, capSteps] = decodeTopologyCandidate(x, mpc, nOpen, cfg);
        txt = sprintf('open=%s caps=%s', mat2str(openLines), mat2str(capSteps));
    catch
        txt = 'open/caps decoding failed';
    end
end

function txt = droopProgressDecode(x, nPV)
    try
        Vsp = x(1:nPV);
        Vband = x(nPV+1:end);
        txt = sprintf('Vsp=%s Vband=%s', mat2str(Vsp,4), mat2str(Vband,4));
    catch
        txt = 'droop decoding failed';
    end
end

function x = reflectBounds(x, lb, ub)
    for j = 1:numel(x)
        lo = lb(j); hi = ub(j);
        if hi <= lo
            x(j) = lo; continue;
        end
        d = hi - lo;
        x(j) = lo + abs(mod(x(j) - hi, 2*d) - d);
        x(j) = min(hi, max(lo, x(j)));
    end
end

function [openLines, capSteps] = decodeTopologyCandidate(x, mpc, nOpen, cfg)
    nbr = size(mpc.branch,1);
    nCaps = numel(cfg.capBuses);

    candidate = round(x(1:nOpen));
    candidate = min(nbr, max(1, candidate));
    openLines = unique(candidate, 'stable');

    % If duplicates or infeasible topology occur, replace with a valid radial set.
    if numel(openLines) ~= nOpen || ~isRadialAfterOpening(mpc, openLines)
        openLines = randomRadialOpenSet(mpc, nOpen);
    end

    capRaw = round(x(nOpen+1:nOpen+nCaps));
    capSteps = min(cfg.capMaxSteps, max(0, capRaw));
end

function openLines = randomRadialOpenSet(mpc, nOpen)
    nbr = size(mpc.branch,1);
    while true
        openLines = sampleK(nbr, nOpen);
        if isRadialAfterOpening(mpc, openLines)
            openLines = sort(openLines);
            return;
        end
    end
end

function idx = sampleK(n, k)
    arr = 1:n;
    idx = zeros(1,k);
    for i = 1:k
        j = randi([i n]);
        tmp = arr(i); arr(i) = arr(j); arr(j) = tmp;
        idx(i) = arr(i);
    end
end

function ok = isRadialAfterOpening(mpc, openLines)
    nb = size(mpc.bus,1);
    nbr = size(mpc.branch,1);
    status = ones(nbr,1);
    status(openLines) = 0;
    B = [mpc.branch(:,1:2), status];
    ok = isRadial(B, nb);
end

function ok = isRadial(branchStatus, nb)
    closed = branchStatus(branchStatus(:,3)==1,1:2);
    if size(closed,1) ~= nb - 1
        ok = false; return;
    end

    adj = cell(nb,1);
    for i = 1:size(closed,1)
        u = closed(i,1); v = closed(i,2);
        if u < 1 || u > nb || v < 1 || v > nb
            ok = false; return;
        end
        adj{u}(end+1) = v; %#ok<AGROW>
        adj{v}(end+1) = u; %#ok<AGROW>
    end

    visited = false(nb,1);
    queue = 1; visited(1) = true; head = 1;
    while head <= numel(queue)
        node = queue(head); head = head + 1;
        for nbrNode = adj{node}
            if ~visited(nbrNode)
                visited(nbrNode) = true;
                queue(end+1) = nbrNode; %#ok<AGROW>
            end
        end
    end
    ok = all(visited);
end

function metrics = computeMetrics(V, cfg)
    under = V < cfg.Vmin;
    over  = V > cfg.Vmax;
    metrics.minV = min(V(:));
    metrics.maxV = max(V(:));
    metrics.meanV = mean(V(:));
    metrics.totalDeviation = sum((V(:) - 1).^2);
    metrics.totalPenalty = sum(max(cfg.Vmin - V(:),0).^2 + max(V(:) - cfg.Vmax,0).^2);
    metrics.underVoltageCount = nnz(under);
    metrics.overVoltageCount = nnz(over);
    metrics.violationCount = nnz(under | over);
end

function exportControlSchedule(controlLog, cfg, fpath)
    T = numel(controlLog);
    TimeStep = (1:T)';
    Mode = strings(T,1);
    OpenLines = strings(T,1);
    CapSteps = strings(T,1);
    Vsp = strings(T,1);
    Vband = strings(T,1);
    Score = nan(T,1);
    for t = 1:T
        Mode(t) = string(controlLog(t).mode);
        OpenLines(t) = string(mat2str(controlLog(t).openLines));
        CapSteps(t) = string(mat2str(controlLog(t).capSteps));
        Vsp(t) = string(mat2str(controlLog(t).Vsp,4));
        Vband(t) = string(mat2str(controlLog(t).Vband,4));
        Score(t) = controlLog(t).score;
    end
    TBL = table(TimeStep, Mode, OpenLines, CapSteps, Vsp, Vband, Score);
    writetable(TBL, fpath);
end

function writeLatexSummaryTable(summaryTable, fpath)
    fid = fopen(fpath, 'w');
    if fid < 0
        warning('Could not write LaTeX summary table.'); return;
    end
    fprintf(fid, '%% Auto-generated by hybrid_dnr_cap_droop_33bus_final_verbose.m\n');
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Voltage-quality comparison for the IEEE 33-bus system.}\n');
    fprintf(fid, '\\label{tab:summary_33bus}\n');
    fprintf(fid, '\\begin{tabular}{lrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Case & Min V & Max V & $\\sum(V-1)^2$ & Viol. count \\\\ \n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(summaryTable)
        fprintf(fid, '%s & %.4f & %.4f & %.4g & %d \\\\ \n', ...
            summaryTable.Case{i}, summaryTable.MinV(i), summaryTable.MaxV(i), ...
            summaryTable.TotalVoltageDeviation(i), summaryTable.TotalViolationCount(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end

function plotAndSaveFigures(cfg, V_no, V_hybrid, controlled, allDroopCurves, pvMW, loadMult)
    t = 1:cfg.T;

    fig = figure('Color','w');
    plot(t, min(V_no,[],2), '--', 'LineWidth', 1.2); hold on;
    plot(t, min(V_hybrid,[],2), '-', 'LineWidth', 1.4);
    yline(cfg.Vmin, ':'); yline(cfg.Vmax, ':');
    grid on; xlabel('Time step'); ylabel('Minimum voltage (p.u.)');
    legend('No control','Hybrid','Location','best');
    saveFigure(fig, fullfile(cfg.outputDir, 'fig_min_voltage_33bus.png'));

    fig = figure('Color','w');
    plot(t, mean(V_no,2), '--', 'LineWidth', 1.2); hold on;
    plot(t, mean(V_hybrid,2), '-', 'LineWidth', 1.4);
    grid on; xlabel('Time step'); ylabel('Mean voltage (p.u.)');
    legend('No control','Hybrid','Location','best');
    saveFigure(fig, fullfile(cfg.outputDir, 'fig_mean_voltage_33bus.png'));

    fig = figure('Color','w');
    imagesc(V_hybrid); colorbar;
    xlabel('Bus'); ylabel('Time step'); title('Hybrid voltage heatmap, IEEE 33-bus');
    saveFigure(fig, fullfile(cfg.outputDir, 'fig_voltage_heatmap_33bus.png'));

    fig = figure('Color','w');
    hold on; grid on;
    for w = 1:numel(controlled.topologyCurve)
        plot(controlled.topologyCurve{w}, 'LineWidth', 1.2);
    end
    xlabel('GWO iteration'); ylabel('Best objective');
    legend(compose('DNR window %d', 1:numel(controlled.topologyCurve)), 'Location','best');
    saveFigure(fig, fullfile(cfg.outputDir, 'fig_gwo_topology_convergence_33bus.png'));

    fig = figure('Color','w');
    plot(t, sum(loadMult,1), 'LineWidth', 1.2); hold on;
    plot(t, sum(pvMW,1), 'LineWidth', 1.2);
    grid on; xlabel('Time step'); ylabel('Aggregate profile');
    legend('Load multiplier sum','PV power (MW)','Location','best');
    saveFigure(fig, fullfile(cfg.outputDir, 'fig_load_pv_profiles_33bus.png'));

    if ~isempty(allDroopCurves)
        fig = figure('Color','w'); hold on; grid on;
        for w = 1:numel(allDroopCurves)
            C = allDroopCurves{w};
            plot(mean(C,1,'omitnan'), 'LineWidth', 1.2);
        end
        xlabel('GWO iteration'); ylabel('Mean best objective across time steps');
        legend(compose('Droop window %d', 1:numel(allDroopCurves)), 'Location','best');
        saveFigure(fig, fullfile(cfg.outputDir, 'fig_gwo_droop_convergence_33bus.png'));
    end
end

function saveFigure(fig, fpath)
    try
        exportgraphics(fig, fpath, 'Resolution', 300);
    catch
        saveas(fig, fpath);
    end
end
