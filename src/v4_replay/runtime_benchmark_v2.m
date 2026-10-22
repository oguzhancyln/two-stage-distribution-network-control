%% runtime_benchmark_v2.m
% Clean runtime benchmark for the IEEE Access manuscript:
% "Condition-Aware Hybrid Voltage Control for PV-Rich Distribution Networks"
%
% Purpose
% -------
% This script measures computational time without repeating the lengthy
% five-run robustness experiment. It separates:
%   1) selected-schedule replay time,
%   2) GWO runtime for scheduled DNR+CAP optimization windows,
%   3) GWO runtime for representative or complete inverter-tuning steps.
%
% IMPORTANT SCIENTIFIC INTERPRETATION
% -----------------------------------
% - Do not claim online feasibility before reviewing the measured values.
% - The DNR+CAP layer may be positioned as offline or intraday supervisory
%   scheduling if its runtime exceeds the 15-minute operating interval.
% - The selected schedule and local inverter actions remain operationally
%   deployable even when schedule generation is performed ahead of time.
%
% Run sequence
% ------------
% A) First use:
%       cfg.mode = 'screening';
%    This measures one DNR+CAP window and two inverter-tuning steps per
%    feeder using the FINAL population and iteration settings.
%
% B) If the screening values are manageable, change:
%       cfg.mode = 'full';
%    This measures all three DNR+CAP windows and all 27 inverter-tuning
%    steps per feeder.
%
% Requirements
% ------------
% - MATPOWER on path
% - case33bw available through MATPOWER
% - case69modified.m or case69.m available on MATLAB path
% - bundled verified profiles and schedules in the same folder
%
% Practical advice
% ----------------
% Use a charger, disable sleep mode, close other heavy applications, and
% run the FULL mode in an uninterrupted session. If Windows restarts or
% the computer sleeps, discard that run and repeat it.

clear; clc; close all;
warning('off','all');

thisFile = mfilename('fullpath');
if isempty(thisFile)
    scriptDir = pwd;
else
    scriptDir = fileparts(thisFile);
end
addpath(scriptDir);
cd(scriptDir);

%% -------------------------- USER SETTINGS ---------------------------
cfg.baseDir = scriptDir;

% FIRST RUN: 'screening'
% PUBLICATION RUNTIME RUN AFTER SCREENING: 'full'
cfg.mode = 'screening';

% One uninterrupted repetition is sufficient for the first assessment.
% Increase to 3 only after the single full run is shown to be manageable.
cfg.nRepeats = 1;

% Final-paper optimizer settings
cfg.topoAgents  = 50;
cfg.topoIter    = 100;
cfg.droopAgents = 40;
cfg.droopIter   = 50;

cfg.verboseIter = false;
cfg.baseSeed = 62001;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambda = 1e5;
cfg.pfmin = 0.90;
cfg.maxCapStep = 5;
cfg.capStepMVAr = 60/1000;  % 60 kVAr per step, expressed in MVAr

cfg.allTopoWindows = [1 24; 40 60; 73 96];
cfg.allDroopSteps = [25:39, 61:72];

switch lower(cfg.mode)
    case 'screening'
        % The late-day DNR+CAP window is measured because it is one of the
        % longest scheduled windows. Two PV-active inverter steps are used.
        cfg.measuredTopoWindows = [73 96];
        cfg.measuredDroopSteps = [30 65];
    case 'full'
        cfg.measuredTopoWindows = cfg.allTopoWindows;
        cfg.measuredDroopSteps = cfg.allDroopSteps;
    otherwise
        error('Unknown cfg.mode: %s. Use screening or full.', cfg.mode);
end

cfg.outputDir = fullfile(cfg.baseDir, sprintf('results_runtime_%s', lower(cfg.mode)));
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end

fprintf('\n============================================================\n');
fprintf(' CLEAN RUNTIME BENCHMARK | mode=%s\n', upper(cfg.mode));
fprintf('============================================================\n');
fprintf('Final settings: DNR+CAP=%d agents x %d iterations | DROOP=%d agents x %d iterations\n', ...
    cfg.topoAgents, cfg.topoIter, cfg.droopAgents, cfg.droopIter);
fprintf('Measured DNR+CAP windows: %s\n', mat2str(cfg.measuredTopoWindows));
fprintf('Measured inverter steps : %s\n', mat2str(cfg.measuredDroopSteps));

writeHardwareInfo(cfg);

systems = makeSystemConfigs(cfg);
detailRows = {};
summaryRows = {};

for sidx = 1:numel(systems)
    sys = systems(sidx);
    fprintf('\n================ %s runtime benchmark ================\n', sys.name);

    [mpc0, Lprof, PVprof, defaultOpen] = loadSystemAndProfiles(sys, cfg);
    sched = readSchedule(sys.scheduleFile, size(Lprof,2));

    fprintf('Loaded %s: %d buses, %d branches, default open=%s\n', ...
        sys.caseName, size(mpc0.bus,1), size(mpc0.branch,1), mat2str(defaultOpen));
    fprintf('Profiles: load range=[%.3f, %.3f], PV max=%.4f MW/site\n', ...
        min(Lprof(:)), max(Lprof(:)), max(PVprof(:)));

    % Warm-up calls prevent one-time MATLAB/MATPOWER initialization from
    % contaminating the reported measurements.
    fprintf('Warm-up power-flow calls ...\n');
    runOneTimeStepNoDroop(mpc0,Lprof,PVprof,sys,cfg,1,defaultOpen,zeros(1,numel(sys.capBuses)));
    warmStep = cfg.measuredDroopSteps(1);
    warmRow = sched(warmStep);
    if warmRow.HasDroop
        etaWarm = [warmRow.Vsp warmRow.Vband];
        runOneTimeStepWithDroop(mpc0,Lprof,PVprof,sys,cfg,warmStep,warmRow.OpenLines,warmRow.CapSteps,etaWarm);
    end

    for rep = 1:cfg.nRepeats
        seed = cfg.baseSeed + 1000*sidx + rep;
        rng(seed,'twister');
        fprintf('\n--- %s repetition %d/%d | seed=%d ---\n',sys.name,rep,cfg.nRepeats,seed);

        % A. Selected-schedule replay timing
        replayTimer = tic;
        [acceptedReplay,rejectedReplay] = replaySelectedSchedule(mpc0,Lprof,PVprof,sys,cfg,sched,defaultOpen);
        replaySec = toc(replayTimer);
        detailRows(end+1,:) = {sys.name,cfg.mode,rep,'SELECTED_SCHEDULE_REPLAY','1..96', ...
            0,0,replaySec,NaN,'','',acceptedReplay,rejectedReplay}; %#ok<SAGROW>
        fprintf('Replay 1..96 DONE | %.3f s | droop accepted=%d rejected=%d\n', ...
            replaySec,acceptedReplay,rejectedReplay);

        % B. Scheduled DNR+CAP optimization timing
        topoTimes = zeros(size(cfg.measuredTopoWindows,1),1);
        for w = 1:size(cfg.measuredTopoWindows,1)
            tRange = cfg.measuredTopoWindows(w,1):cfg.measuredTopoWindows(w,2);
            label = sprintf('%s runtime TOPO+CAP t=%d..%d',sys.name,tRange(1),tRange(end));
            taskTimer = tic;
            [bestOpen,bestCaps,bestObj] = gwoTopologyCap(mpc0,Lprof,PVprof,sys,cfg,tRange,label);
            sec = toc(taskTimer);
            topoTimes(w) = sec;
            detailRows(end+1,:) = {sys.name,cfg.mode,rep,'DNR_CAP_OPT',sprintf('%d..%d',tRange(1),tRange(end)), ...
                cfg.topoAgents,cfg.topoIter,sec,bestObj,mat2str(bestOpen),mat2str(bestCaps),NaN,NaN}; %#ok<SAGROW>
            fprintf('DNR+CAP %d..%d DONE | %.1f s | best=%.6g | open=%s caps=%s\n', ...
                tRange(1),tRange(end),sec,bestObj,mat2str(bestOpen),mat2str(bestCaps));
        end

        % C. Inverter-tuning timing using retained selected-schedule states
        droopTimes = zeros(numel(cfg.measuredDroopSteps),1);
        for d = 1:numel(cfg.measuredDroopSteps)
            t = cfg.measuredDroopSteps(d);
            row = sched(t);
            if isempty(row.OpenLines), row.OpenLines = defaultOpen; end
            if isempty(row.CapSteps), row.CapSteps = zeros(1,numel(sys.capBuses)); end

            label = sprintf('%s runtime DROOP t=%d',sys.name,t);
            taskTimer = tic;
            [~,bestObj,~,~] = gwoDroopOneStep(mpc0,Lprof,PVprof,sys,cfg,t,row.OpenLines,row.CapSteps,label);
            sec = toc(taskTimer);
            droopTimes(d) = sec;
            detailRows(end+1,:) = {sys.name,cfg.mode,rep,'INVERTER_TUNING_OPT',sprintf('%d',t), ...
                cfg.droopAgents,cfg.droopIter,sec,bestObj,mat2str(row.OpenLines),mat2str(row.CapSteps),NaN,NaN}; %#ok<SAGROW>
            fprintf('Inverter tuning t=%d DONE | %.1f s | best=%.6g\n',t,sec,bestObj);
        end

        totalMeasured = replaySec + sum(topoTimes) + sum(droopTimes);
        summaryRows(end+1,:) = {sys.name,cfg.mode,rep,seed,replaySec, ...
            sum(topoTimes),mean(topoTimes),max(topoTimes), ...
            sum(droopTimes),mean(droopTimes),max(droopTimes),numel(droopTimes), ...
            totalMeasured,max(topoTimes)<=900,max(droopTimes)<=900, ...
            acceptedReplay,rejectedReplay}; %#ok<SAGROW>

        fprintf('%s repetition %d SUMMARY | replay=%.2fs | DNR+CAP total=%.1fs mean=%.1fs max=%.1fs | ', ...
            sys.name,rep,replaySec,sum(topoTimes),mean(topoTimes),max(topoTimes));
        fprintf('inverter total=%.1fs mean=%.1fs max=%.1fs | measured total=%.1fs\n', ...
            sum(droopTimes),mean(droopTimes),max(droopTimes),totalMeasured);
    end
end

detailNames = {'System','Mode','Repetition','Task','WindowOrStep','Agents','Iterations', ...
    'RuntimeSeconds','BestObjective','OpenLines','CapSteps','ReplayAccepted','ReplayRejected'};
summaryNames = {'System','Mode','Repetition','Seed','Replay96Steps_s', ...
    'DNRCAP_Total_s','DNRCAP_MeanWindow_s','DNRCAP_MaxWindow_s', ...
    'Inverter_Total_s','Inverter_MeanStep_s','Inverter_MaxStep_s','MeasuredInverterSteps', ...
    'TotalMeasured_s','DNRCAP_MaxUnder15min','Inverter_MaxUnder15min', ...
    'ReplayAccepted','ReplayRejected'};

Tdetail = cell2table(detailRows,'VariableNames',detailNames);
Tsummary = cell2table(summaryRows,'VariableNames',summaryNames);

writetable(Tdetail,fullfile(cfg.outputDir,'runtime_detail.csv'));
writetable(Tsummary,fullfile(cfg.outputDir,'runtime_summary.csv'));
writeLatexRuntimeTable(Tsummary,fullfile(cfg.outputDir,'runtime_summary_latex_table.tex'));
writeInterpretationTemplate(Tsummary,cfg,fullfile(cfg.outputDir,'runtime_interpretation_template.txt'));

fprintf('\nDone. Runtime outputs saved under:\n  %s\n',cfg.outputDir);
fprintf('\nNEXT ACTION:\n');
if strcmpi(cfg.mode,'screening')
    fprintf('Send runtime_summary.csv. Do NOT switch to full mode until the screening times are reviewed.\n');
else
    fprintf('Send runtime_summary.csv and hardware_info.txt for manuscript integration.\n');
end

%% ============================ LOCAL FUNCTIONS ============================
function systems = makeSystemConfigs(cfg)
    systems(1).name = '33bus';
    systems(1).caseName = 'case33bw';
    systems(1).fallbackCaseName = 'case33bw';
    systems(1).defaultOpen = [33 34 35 36 37];
    systems(1).pvBuses = [21 24 32 17 18];
    systems(1).capBuses = [18 25 33];
    systems(1).pvRatedMW = 0.20;
    systems(1).loadFile = locateInputFile(cfg.baseDir,'myloadspvs.mat');
    systems(1).pvFile = locateInputFile(cfg.baseDir,'myloadspvs.mat');
    systems(1).scheduleFile = locateInputFile(cfg.baseDir,'control_schedule_33bus.csv');

    systems(2).name = '69bus';
    systems(2).caseName = 'case69modified';
    systems(2).fallbackCaseName = 'case69';
    systems(2).defaultOpen = [69 70 71 72 73];
    systems(2).pvBuses = [15 18 27 35 46 60];
    systems(2).capBuses = [7 23 30 39 55];
    systems(2).pvRatedMW = 1.00;
    systems(2).loadFile = locateInputFile(cfg.baseDir,'myloads.mat');
    systems(2).pvFile = locateInputFile(cfg.baseDir,'mypvs.mat');
    systems(2).scheduleFile = locateInputFile(cfg.baseDir,'control_schedule_69bus.csv');
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
        error('%s: required verified profile files could not be loaded. Do not use synthetic profiles for the runtime benchmark.', sys.name);
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


function filepath = locateInputFile(baseDir,filename)
    % Resolve bundled publication inputs without depending on MATLAB's
    % current working directory. Direct script-folder matches are preferred.
    directCandidates = {
        fullfile(baseDir,filename)
        fullfile(pwd,filename)
    };

    for k = 1:numel(directCandidates)
        if exist(directCandidates{k},'file')
            filepath = directCandidates{k};
            fprintf('Input resolved: %s\n',filepath);
            return;
        end
    end

    % Search recursively below the script folder. This supports extracting
    % the ZIP into a subfolder while launching the script from its parent.
    recursiveMatches = dir(fullfile(baseDir,'**',filename));
    recursiveMatches = recursiveMatches(~[recursiveMatches.isdir]);

    if isempty(recursiveMatches)
        error(['Required input file not found: %s\n' ...
               'Place it in the same folder as runtime_benchmark_v2.m, ' ...
               'or extract the complete ZIP bundle into one folder.'], filename);
    end

    fullMatches = arrayfun(@(x) fullfile(x.folder,x.name), recursiveMatches, ...
        'UniformOutput',false);
    fullMatches = unique(fullMatches);

    if numel(fullMatches) > 1
        fprintf('\nMultiple copies of %s were found:\n',filename);
        for j = 1:numel(fullMatches)
            fprintf('  %s\n',fullMatches{j});
        end
        error(['Multiple copies of %s were found. Keep one verified copy ' ...
               'in the script folder to avoid using stale inputs.'],filename);
    end

    filepath = fullMatches{1};
    fprintf('Input resolved recursively: %s\n',filepath);
end

function sched = readSchedule(filename,T)
    if ~exist(filename,'file')
        error('Schedule file not found: %s',filename);
    end
    TBL = readtable(filename,'TextType','string');
    sched = repmat(struct('OpenLines',[],'CapSteps',[],'Vsp',[],'Vband',[],'HasDroop',false),T,1);
    for i = 1:height(TBL)
        t = TBL.TimeStep(i);
        if t<1 || t>T, continue; end
        sched(t).OpenLines = parseVector(TBL.OpenLines(i));
        sched(t).CapSteps  = parseVector(TBL.CapSteps(i));
        sched(t).Vsp       = parseVector(TBL.Vsp(i));
        sched(t).Vband     = parseVector(TBL.Vband(i));
        sched(t).HasDroop  = ~isempty(sched(t).Vsp) && ~isempty(sched(t).Vband);
    end
end

function v = parseVector(x)
    if ismissing(x), v=[]; return; end
    s = strtrim(char(x));
    if isempty(s) || strcmp(s,'[]'), v=[]; return; end
    s = strrep(s,'[',''); s = strrep(s,']',''); s = strrep(s,',',' ');
    v = sscanf(s,'%f').';
end

function [accepted,rejected] = replaySelectedSchedule(mpc0,Lprof,PVprof,sys,cfg,sched,defaultOpen)
    T = size(Lprof,2);
    accepted = 0; rejected = 0;
    for t = 1:T
        row = sched(t);
        if isempty(row.OpenLines), row.OpenLines = defaultOpen; end
        if isempty(row.CapSteps), row.CapSteps = zeros(1,numel(sys.capBuses)); end
        if row.HasDroop
            eta = [row.Vsp row.Vband];
            [~,~,Jpre,Jpost] = runOneTimeStepWithDroop(mpc0,Lprof,PVprof,sys,cfg,t,row.OpenLines,row.CapSteps,eta);
            if Jpost <= Jpre
                accepted = accepted + 1;
            else
                rejected = rejected + 1;
            end
        else
            runOneTimeStepNoDroop(mpc0,Lprof,PVprof,sys,cfg,t,row.OpenLines,row.CapSteps);
        end
    end
end

function writeHardwareInfo(cfg)
    outfile = fullfile(cfg.outputDir,'hardware_info.txt');
    fid = fopen(outfile,'w');
    fprintf(fid,'Clean runtime benchmark hardware information\n');
    fprintf(fid,'==========================================\n');
    fprintf(fid,'MATLAB version: %s\n',version);
    fprintf(fid,'Computer identifier: %s\n',computer);
    fprintf(fid,'Benchmark mode: %s\n',cfg.mode);
    fprintf(fid,'DNR+CAP settings: %d agents x %d iterations\n',cfg.topoAgents,cfg.topoIter);
    fprintf(fid,'Inverter-tuning settings: %d agents x %d iterations\n',cfg.droopAgents,cfg.droopIter);
    try
        fprintf(fid,'MATLAB detected cores: %d\n',feature('numcores'));
    catch
        fprintf(fid,'MATLAB detected cores: unavailable\n');
    end
    try
        [status,cpu] = system('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor).Name"');
        if status==0, fprintf(fid,'CPU: %s\n',strtrim(cpu)); end
    catch
        fprintf(fid,'CPU: unavailable\n');
    end
    try
        [status,mem] = system('powershell -NoProfile -Command "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)"');
        if status==0, fprintf(fid,'Installed RAM (GB): %s\n',strtrim(mem)); end
    catch
        fprintf(fid,'Installed RAM: unavailable\n');
    end
    fprintf(fid,'\nManual note: discard this run if the computer slept, restarted, or was heavily used during the benchmark.\n');
    fclose(fid);
end

function writeLatexRuntimeTable(T,outFile)
    fid = fopen(outFile,'w');
    fprintf(fid,'%% Auto-generated by runtime_benchmark_v2.m\n');
    fprintf(fid,'\\begin{table*}[!t]\\n\\centering\\n');
    fprintf(fid,'\\caption{Clean stage-decomposed runtime benchmark under the final optimizer settings.}\\n');
    fprintf(fid,'\\label{tab:runtime_benchmark}\\n\\scriptsize\\n');
    fprintf(fid,'\\begin{tabular}{lrrrrrr}\\n\\hline\\n');
    fprintf(fid,'System & Replay 96 steps (s) & DNR+CAP total (s) & DNR+CAP max window (s) & Inverter mean step (s) & Inverter max step (s) & Measured total (s) \\\\\\\\ \\hline\\n');
    for i=1:height(T)
        fprintf(fid,'%s & %.2f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\\\\\ \\n', ...
            T.System{i},T.Replay96Steps_s(i),T.DNRCAP_Total_s(i),T.DNRCAP_MaxWindow_s(i), ...
            T.Inverter_MeanStep_s(i),T.Inverter_MaxStep_s(i),T.TotalMeasured_s(i));
    end
    fprintf(fid,'\\hline\\n\\end{tabular}\\n\\end{table*}\\n');
    fclose(fid);
end

function writeInterpretationTemplate(T,cfg,outFile)
    fid = fopen(outFile,'w');
    fprintf(fid,'RUNTIME BENCHMARK INTERPRETATION TEMPLATE\\n');
    fprintf(fid,'========================================\\n\\n');
    fprintf(fid,'Mode: %s\\n\\n',cfg.mode);
    fprintf(fid,'Do not copy screening values into the paper. Use FULL mode for manuscript reporting.\\n\\n');
    for i=1:height(T)
        fprintf(fid,'%s\\n',T.System{i});
        fprintf(fid,'- selected-schedule replay: %.2f s\\n',T.Replay96Steps_s(i));
        fprintf(fid,'- DNR+CAP maximum measured window: %.1f s\\n',T.DNRCAP_MaxWindow_s(i));
        fprintf(fid,'- inverter-tuning mean measured step: %.1f s\\n',T.Inverter_MeanStep_s(i));
        fprintf(fid,'- inverter-tuning maximum measured step: %.1f s\\n',T.Inverter_MaxStep_s(i));
        fprintf(fid,'- DNR+CAP maximum below 15 minutes: %d\\n',T.DNRCAP_MaxUnder15min(i));
        fprintf(fid,'- inverter-tuning maximum below 15 minutes: %d\\n\\n',T.Inverter_MaxUnder15min(i));
    end
    fprintf(fid,'Interpretation rule:\\n');
    fprintf(fid,'- If DNR+CAP exceeds 900 s, describe it as offline or intraday supervisory schedule generation.\\n');
    fprintf(fid,'- Do not claim real-time optimization solely because the replay is fast.\\n');
    fprintf(fid,'- Report hardware information together with FULL-mode timings.\\n');
    fclose(fid);
end
