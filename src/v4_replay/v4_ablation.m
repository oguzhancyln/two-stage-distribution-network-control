%% run_numeric_ablation_33_69_from_schedules_v4_consistent_qmax.m
% Numeric ablation generator for the IEEE 33-bus and 69-bus paper cases.
% V4 fixes the replay inverter-Q capability model so that it matches the
% manuscript and the full multi-run workflow.
%
% It evaluates four cases using the final schedules already produced by the
% safe hybrid scripts:
%   1) No control
%   2) DNR+CAP only
%   3) DROOP only
%   4) Proposed hybrid DNR+CAP+DROOP
%
% Fast default behavior:
%   - DNR+CAP only uses the final open-line and capacitor schedules.
%   - DROOP only uses the final per-time-step Vsp/Vband droop settings from
%     the schedule but applies them on the default radial topology without
%     capacitor switching. This is fast and reproducible.
%   - Proposed hybrid replays the final schedule with a droop safety gate.
%     The safety gate accepts inverter droop only if it improves the
%     instantaneous voltage objective relative to the pre-droop state. This
%     prevents droop replay from overcompensating under schedule/profile
%     reconstruction mismatches and matches the safe operational logic used
%     in the final runs.
%
% Optional fairer but slower behavior:
%   Set cfg.reoptimizeDroopOnly = true below. Then DROOP only will optimize
%   Vsp/Vband again under the default radial topology and zero capacitors for
%   each scheduled droop time step. This is more expensive but stronger for
%   publication.
%
% Dependencies:
%   - MATPOWER on MATLAB path
%   - case33bw and case69modified available to loadcase()
%   - control_schedule_33bus.csv, control_schedule_69bus.csv
%   - myloadspvs.mat for 33-bus, with myloads and myPVs
%   - myloads.mat and mypvs.mat for 69-bus, with myloads and PV
%
% Output:
%   results_ablation_33_69_v4_consistent_qmax/
%      ablation_summary_33bus.csv
%      ablation_summary_69bus.csv
%      ablation_summary_combined.csv
%      ablation_summary_33bus_latex_table.tex
%      ablation_summary_69bus_latex_table.tex
%      ablation_results_33bus.mat
%      ablation_results_69bus.mat
%
% Consistent-Qmax ablation helper, June 2026.

clear; clc; close all;
warning('off','all');
rng(42);

% Resolve bundled input files relative to this script rather than the
% current MATLAB working directory. This prevents missing-file errors when
% the script is launched from another folder.
scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end
addpath(scriptDir);

%% ---------------------------- USER SETTINGS ----------------------------
cfg.baseDir = scriptDir;
cfg.outputDir = fullfile(cfg.baseDir, 'results_ablation_33_69_v4_consistent_qmax');
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end

cfg.T       = 96;
cfg.dtHours = 0.25;
cfg.Vmin    = 0.95;
cfg.Vmax    = 1.05;
cfg.lambdaUnder = 1e5;
cfg.lambdaOver  = 1e5;
cfg.PFmin        = 0.90;

% Safe operational replay: apply droop only if it improves the current
% time-step objective compared with the pre-droop condition.
% Keep true for publication-grade ablation; otherwise a stored droop schedule
% can overcompensate when replayed under a slightly different profile/scaling.
cfg.useDroopSafetyGate = true;

% Fast/reproducible default. For a stronger but slower paper run, set true.
cfg.reoptimizeDroopOnly = true;

% If cfg.reoptimizeDroopOnly=true, these settings are used.
cfg.gwoDroopAgents = 30;
cfg.gwoDroopIter   = 30;
cfg.verboseGWO     = true;

% Droop parameter bounds.
cfg.VspLB   = 0.95;
cfg.VspUB   = 1.05;
cfg.VbandLB = 0.01;
cfg.VbandUB = 0.05;

%% --------------------------- RUN BOTH CASES ----------------------------
case33 = struct();
case33.name          = '33bus';
case33.caseName      = 'case33bw';
case33.scheduleFile  = fullfile(cfg.baseDir, 'control_schedule_33bus.csv');
case33.profileType   = 'combined33';       % myloadspvs.mat: myloads, myPVs
case33.profileFile   = fullfile(cfg.baseDir, 'myloadspvs.mat');
case33.pvBuses       = [21 24 32 17 18];
case33.capBuses      = [18 25 33];
case33.capStepMVAr   = 60/1000;
case33.pvMode        = 'normalized_to_rating';
case33.pvRatedMW     = 0.20;              % must match final 33-bus run

case69 = struct();
case69.name          = '69bus';
case69.caseName      = 'case69modified';
case69.scheduleFile  = fullfile(cfg.baseDir, 'control_schedule_69bus.csv');
case69.profileType   = 'split69';          % myloads.mat + mypvs.mat: PV
case69.loadFile      = fullfile(cfg.baseDir, 'myloads.mat');
case69.pvFile        = fullfile(cfg.baseDir, 'mypvs.mat');
case69.pvBuses       = [15 18 27 35 46 60];
case69.capBuses      = [7 23 30 39 55];
case69.capStepMVAr   = 60/1000;
case69.pvMode        = 'MW';
case69.pvRatedMW     = 1.00;              % verified peak MW per PV site

fprintf('\n================ IEEE 33-bus numeric ablation ================\n');
res33 = runOneSystemAblation(case33, cfg);

fprintf('\n================ IEEE 69-bus numeric ablation ================\n');
res69 = runOneSystemAblation(case69, cfg);

% Combined summary
combined = [res33.summary; res69.summary];
writetable(combined, fullfile(cfg.outputDir,'ablation_summary_combined.csv'));

fprintf('\nDone. Numeric ablation files saved under: %s\n', cfg.outputDir);

%% =======================================================================
%% Local functions
%% =======================================================================
function result = runOneSystemAblation(sys, cfg)
    mpc0 = loadcase(sys.caseName);
    nb   = size(mpc0.bus,1);
    nbr  = size(mpc0.branch,1);
    nOpen = nbr - (nb - 1);
    if nOpen <= 0
        error('%s: case has no reconfiguration/tie lines. nOpen=%d', sys.name, nOpen);
    end

    baseOpen = determineDefaultOpenLines(mpc0, nOpen);
    fprintf('Loaded %s: %d buses, %d branches, default open=%s\n', sys.caseName, nb, nbr, mat2str(baseOpen));

    [loadMult, Ppv] = loadProfiles(sys, cfg.T, nb, numel(sys.pvBuses));
    fprintf('Profiles loaded: load range=[%.3f, %.3f], PV max=%.4f MW/site\n', ...
        min(loadMult(:)), max(loadMult(:)), max(Ppv(:)));

    sched = readSchedule(sys.scheduleFile, cfg.T);

    % Evaluate four ablation cases.
    cases = {'No control','DNR+CAP only','DROOP only','Proposed hybrid'};
    Vcase = struct();
    scoreCase = struct();

    for c = 1:numel(cases)
        caseName = cases{c};
        fprintf('\nEvaluating %s / %s ...\n', sys.name, caseName);
        Vmat = nan(cfg.T, nb);
        scoreVec = nan(cfg.T,1);

        for t = 1:cfg.T
            row = sched(t);
            switch caseName
                case 'No control'
                    openLines = baseOpen;
                    capSteps  = zeros(1, numel(sys.capBuses));
                    useDroop  = false;
                    Vsp = []; Vband = [];

                case 'DNR+CAP only'
                    openLines = row.OpenLines;
                    if isempty(openLines), openLines = baseOpen; end
                    capSteps  = row.CapSteps;
                    if isempty(capSteps), capSteps = zeros(1,numel(sys.capBuses)); end
                    useDroop  = false;
                    Vsp = []; Vband = [];

                case 'DROOP only'
                    openLines = baseOpen;
                    capSteps  = zeros(1, numel(sys.capBuses));
                    if cfg.reoptimizeDroopOnly && row.HasDroop
                        [Vsp, Vband] = optimizeDroopForTime(mpc0, openLines, capSteps, t, loadMult, Ppv, sys, cfg, ...
                            sprintf('%s DROOP-only t=%d', sys.name, t));
                        useDroop = true;
                    else
                        Vsp = row.Vsp; Vband = row.Vband;
                        useDroop = row.HasDroop;
                    end

                case 'Proposed hybrid'
                    openLines = row.OpenLines;
                    if isempty(openLines), openLines = baseOpen; end
                    capSteps  = row.CapSteps;
                    if isempty(capSteps), capSteps = zeros(1,numel(sys.capBuses)); end
                    Vsp       = row.Vsp;
                    Vband     = row.Vband;
                    useDroop  = row.HasDroop;
            end

            [V, score] = simulateTimeStep(mpc0, openLines, capSteps, useDroop, Vsp, Vband, ...
                t, loadMult, Ppv, sys, cfg);
            Vmat(t,:) = V;
            scoreVec(t) = score;
        end

        key = matlab.lang.makeValidName(caseName);
        Vcase.(key) = Vmat;
        scoreCase.(key) = scoreVec;
    end

    % Create summary table.
    S = table();
    S.System = repmat({sys.name}, numel(cases), 1);
    S.Case = cases(:);
    for c = 1:numel(cases)
        key = matlab.lang.makeValidName(cases{c});
        metrics = computeMetrics(Vcase.(key), cfg);
        S.MinV(c,1) = metrics.MinV;
        S.MaxV(c,1) = metrics.MaxV;
        S.MeanV(c,1) = metrics.MeanV;
        S.TotalVoltageDeviation(c,1) = metrics.TotalVoltageDeviation;
        S.TotalViolationPenalty(c,1) = metrics.TotalViolationPenalty;
        S.UnderVoltageCount(c,1) = metrics.UnderVoltageCount;
        S.OverVoltageCount(c,1) = metrics.OverVoltageCount;
        S.TotalViolationCount(c,1) = metrics.TotalViolationCount;
    end

    % Improvements relative to no-control.
    baseJ = S.TotalVoltageDeviation(1);
    baseViol = S.TotalViolationCount(1);
    S.DevReductionPercent = 100*(baseJ - S.TotalVoltageDeviation)./max(baseJ, eps);
    S.ViolationReductionPercent = 100*(baseViol - S.TotalViolationCount)./max(baseViol, eps);

    % Save outputs.
    outCsv = fullfile(cfg.outputDir, sprintf('ablation_summary_%s.csv', sys.name));
    writetable(S, outCsv);
    writeLatexAblationTable(S, fullfile(cfg.outputDir, sprintf('ablation_summary_%s_latex_table.tex', sys.name)), sys.name);

    result = struct();
    result.summary = S;
    result.Vcase = Vcase;
    result.scoreCase = scoreCase;
    result.schedule = sched;
    result.baseOpenLines = baseOpen;
    save(fullfile(cfg.outputDir, sprintf('ablation_results_%s.mat', sys.name)), 'result', 'sys', 'cfg');

    disp(S);
end

function [loadMult, Ppv] = loadProfiles(sys, T, nb, nPV)
    switch sys.profileType
        case 'combined33'
            if ~exist(sys.profileFile,'file')
                error('Profile file not found: %s', sys.profileFile);
            end
            data = load(sys.profileFile);
            if ~isfield(data,'myloads') || ~isfield(data,'myPVs')
                error('%s must contain myloads and myPVs.', sys.profileFile);
            end
            loadMult = data.myloads;
            Ppv = data.myPVs;

        case 'split69'
            if ~exist(sys.loadFile,'file') || ~exist(sys.pvFile,'file')
                error('Profile files not found: %s / %s', sys.loadFile, sys.pvFile);
            end
            L = load(sys.loadFile);
            P = load(sys.pvFile);
            if ~isfield(L,'myloads')
                error('%s must contain myloads.', sys.loadFile);
            end
            if ~isfield(P,'PV')
                error('%s must contain PV.', sys.pvFile);
            end
            loadMult = L.myloads;
            Ppv = P.PV;
        otherwise
            error('Unknown profileType: %s', sys.profileType);
    end

    if size(loadMult,1) ~= nb || size(loadMult,2) < T
        error('Load profile size mismatch. Expected %d x %d, got %d x %d.', nb, T, size(loadMult,1), size(loadMult,2));
    end
    if size(Ppv,1) ~= nPV || size(Ppv,2) < T
        error('PV profile size mismatch. Expected %d x %d, got %d x %d.', nPV, T, size(Ppv,1), size(Ppv,2));
    end
    loadMult = loadMult(:,1:T);
    Ppv = Ppv(:,1:T);

    % Ensure profile scaling is consistent with the final optimization runs.
    % For the 33-bus case, myloadspvs.mat stores PV as a normalized/profile
    % matrix in some runs; the final accepted 33-bus optimization used
    % 0.20 MW peak per PV site. Therefore we rescale the profile here so
    % that max(Ppv)=0.20 MW/site.
    if isfield(sys,'pvMode') && strcmpi(sys.pvMode,'normalized_to_rating')
        mx = max(Ppv(:));
        if mx <= 0
            error('PV profile has non-positive maximum; cannot rescale.');
        end
        Ppv = sys.pvRatedMW * Ppv / mx;
        fprintf('PV profile rescaled to %.3f MW/site peak for %s.\n', sys.pvRatedMW, sys.name);
    end
end

function sched = readSchedule(filename, T)
    if ~exist(filename,'file')
        error('Schedule file not found: %s', filename);
    end
    TBL = readtable(filename, 'TextType', 'string');
    sched = repmat(struct('OpenLines',[],'CapSteps',[],'Vsp',[],'Vband',[],'HasDroop',false,'Mode',''), T, 1);
    for i = 1:height(TBL)
        t = TBL.TimeStep(i);
        if t < 1 || t > T, continue; end
        sched(t).Mode = char(TBL.Mode(i));
        sched(t).OpenLines = parseVector(TBL.OpenLines(i));
        sched(t).CapSteps  = parseVector(TBL.CapSteps(i));
        sched(t).Vsp       = parseVector(TBL.Vsp(i));
        sched(t).Vband     = parseVector(TBL.Vband(i));
        sched(t).HasDroop  = ~isempty(sched(t).Vsp) && ~isempty(sched(t).Vband);
    end
end

function v = parseVector(x)
    if ismissing(x)
        v = []; return;
    end
    s = char(x);
    s = strtrim(s);
    if isempty(s) || strcmp(s,'[]')
        v = []; return;
    end
    s = strrep(s,'[','');
    s = strrep(s,']','');
    s = strrep(s,',',' ');
    v = sscanf(s,'%f').';
end

function [Vfinal, score] = simulateTimeStep(mpc0, openLines, capSteps, useDroop, Vsp, Vband, t, loadMult, Ppv, sys, cfg)
    m0 = prepareCase(mpc0, openLines, capSteps, t, loadMult, Ppv, sys);
    mpopt = mpoption('verbose',0,'out.all',0);

    if useDroop
        % First solve the pre-droop state. This is also the fall-back state
        % if the droop action overcompensates or causes a worse objective.
        r0 = runpf(m0, mpopt);
        if ~r0.success
            Vfinal = nan(1,size(mpc0.bus,1)); score = inf; return;
        end
        V0 = abs(r0.bus(:,8)).';
        score0 = objectiveVoltage(V0, cfg);

        % Apply inverter Volt/VAr commands around the measured pre-droop voltage.
        QmaxStep = computeQmax(Ppv(:,t), mpc0.baseMVA, cfg, sys);
        m1 = m0;
        for k = 1:numel(sys.pvBuses)
            idxBus = sys.pvBuses(k);
            Qpu = pvControlSetPointLocal(V0(idxBus), Vsp(k), Vband(k), QmaxStep(k));
            m1.bus(idxBus,4) = m1.bus(idxBus,4) - Qpu*mpc0.baseMVA;
        end
        r1 = runpf(m1, mpopt);
        if ~r1.success
            Vfinal = V0; score = score0; return;
        end
        V1 = abs(r1.bus(:,8)).';
        score1 = objectiveVoltage(V1, cfg);

        % Safety gate: accept droop only if it improves or preserves the
        % instantaneous objective. This avoids harmful overcompensation when
        % replaying a stored droop schedule in ablation experiments.
        if isfield(cfg,'useDroopSafetyGate') && cfg.useDroopSafetyGate && score1 > score0
            Vfinal = V0;
            score = score0;
            return;
        else
            Vfinal = V1;
            score = score1;
            return;
        end
    else
        r = runpf(m0, mpopt);
        if ~r.success
            Vfinal = nan(1,size(mpc0.bus,1)); score = inf; return;
        end
        Vfinal = abs(r.bus(:,8)).';
        score = objectiveVoltage(Vfinal, cfg);
        return;
    end
end

function m = prepareCase(mpc0, openLines, capSteps, t, loadMult, Ppv, sys)
    m = mpc0;
    m.branch(:,11) = 1;
    m.branch(openLines,11) = 0;
    P0 = m.bus(:,3);
    Q0 = m.bus(:,4);
    m.bus(:,3) = P0 .* loadMult(:,t);
    m.bus(:,4) = Q0 .* loadMult(:,t);
    for k = 1:numel(sys.pvBuses)
        m.bus(sys.pvBuses(k),3) = m.bus(sys.pvBuses(k),3) - Ppv(k,t);
    end
    if nargin >= 3 && ~isempty(capSteps)
        for k = 1:min(numel(sys.capBuses), numel(capSteps))
            m.bus(sys.capBuses(k),4) = m.bus(sys.capBuses(k),4) - capSteps(k) * sys.capStepMVAr;
        end
    end
end

function QmaxStep = computeQmax(PactMW, baseMVA, cfg, sys)
    % Consistent dynamic inverter capability model used in the manuscript:
    % Qmax = min(sqrt(S^2-P^2), P*tan(acos(PFmin))).
    if ~isfield(sys,'pvRatedMW') || ~isfinite(sys.pvRatedMW) || sys.pvRatedMW <= 0
        error('%s: sys.pvRatedMW must be a positive finite value.', sys.name);
    end
    SmaxMW = sys.pvRatedMW ./ cfg.PFmin;
    qCircleMVAr = sqrt(max(SmaxMW.^2 - PactMW.^2, 0));
    qPfMVAr = PactMW .* tan(acos(cfg.PFmin));
    qMaxMVAr = min(qCircleMVAr, qPfMVAr);
    QmaxStep = qMaxMVAr ./ baseMVA;
end

function Qcmd = pvControlSetPointLocal(Vmeas, Vsp, Vband, Qmax)
    K = Qmax ./ Vband;
    dV = Vmeas - Vsp;
    Qcmd = -K .* dV;
    if dV <= -Vband
        Qcmd = Qmax;
    elseif dV >= Vband
        Qcmd = -Qmax;
    end
end

function J = objectiveVoltage(V, cfg)
    if any(isnan(V))
        J = inf; return;
    end
    under = max(cfg.Vmin - V, 0);
    over  = max(V - cfg.Vmax, 0);
    J = sum((V - 1).^2) + cfg.lambdaUnder*sum(under.^2) + cfg.lambdaOver*sum(over.^2);
end

function metrics = computeMetrics(Vmat, cfg)
    V = Vmat(:);
    V = V(~isnan(V));
    under = max(cfg.Vmin - V, 0);
    over  = max(V - cfg.Vmax, 0);
    metrics.MinV = min(V);
    metrics.MaxV = max(V);
    metrics.MeanV = mean(V);
    metrics.TotalVoltageDeviation = sum((V - 1).^2);
    metrics.TotalViolationPenalty = sum(under.^2 + over.^2);
    metrics.UnderVoltageCount = sum(V < cfg.Vmin);
    metrics.OverVoltageCount = sum(V > cfg.Vmax);
    metrics.TotalViolationCount = metrics.UnderVoltageCount + metrics.OverVoltageCount;
end

function openLines = determineDefaultOpenLines(mpc, nOpen)
    status = mpc.branch(:,11);
    if numel(status) >= 1
        existingOpen = find(status == 0).';
    else
        existingOpen = [];
    end
    if numel(existingOpen) == nOpen
        openLines = existingOpen;
    else
        openLines = (size(mpc.branch,1)-nOpen+1):size(mpc.branch,1);
    end
end

function [bestVsp, bestVband] = optimizeDroopForTime(mpc0, openLines, capSteps, t, loadMult, Ppv, sys, cfg, label)
    nPV = numel(sys.pvBuses);
    dim = 2*nPV;
    lb = [cfg.VspLB*ones(1,nPV), cfg.VbandLB*ones(1,nPV)];
    ub = [cfg.VspUB*ones(1,nPV), cfg.VbandUB*ones(1,nPV)];
    X = rand(cfg.gwoDroopAgents, dim).*(ub-lb) + lb;
    alphaPos = zeros(1,dim); betaPos = alphaPos; deltaPos = alphaPos;
    alphaScore = inf; betaScore = inf; deltaScore = inf;
    tStart = tic;

    for it = 1:cfg.gwoDroopIter
        feasible = 0;
        for i = 1:cfg.gwoDroopAgents
            X(i,:) = min(max(X(i,:),lb),ub);
            Vsp = X(i,1:nPV);
            Vband = X(i,nPV+1:end);
            [~, score] = simulateTimeStep(mpc0, openLines, capSteps, true, Vsp, Vband, t, loadMult, Ppv, sys, cfg);
            if isfinite(score), feasible = feasible + 1; end
            if score < alphaScore
                deltaScore = betaScore; deltaPos = betaPos;
                betaScore = alphaScore; betaPos = alphaPos;
                alphaScore = score; alphaPos = X(i,:);
            elseif score < betaScore
                deltaScore = betaScore; deltaPos = betaPos;
                betaScore = score; betaPos = X(i,:);
            elseif score < deltaScore
                deltaScore = score; deltaPos = X(i,:);
            end
        end

        a = 2 - it*(2/cfg.gwoDroopIter);
        for i = 1:cfg.gwoDroopAgents
            for j = 1:dim
                r1=rand; r2=rand; A1=2*a*r1-a; C1=2*r2;
                Dalpha=abs(C1*alphaPos(j)-X(i,j)); X1=alphaPos(j)-A1*Dalpha;
                r1=rand; r2=rand; A2=2*a*r1-a; C2=2*r2;
                Dbeta=abs(C2*betaPos(j)-X(i,j)); X2=betaPos(j)-A2*Dbeta;
                r1=rand; r2=rand; A3=2*a*r1-a; C3=2*r2;
                Ddelta=abs(C3*deltaPos(j)-X(i,j)); X3=deltaPos(j)-A3*Ddelta;
                X(i,j)=(X1+X2+X3)/3;
            end
        end
        X = min(max(X,lb),ub);
        if cfg.verboseGWO
            fprintf('[%s] iter %3d/%3d | best=%g | feasible=%d/%d | elapsed=%.1fs\n', ...
                label, it, cfg.gwoDroopIter, alphaScore, feasible, cfg.gwoDroopAgents, toc(tStart));
        end
    end
    bestVsp = alphaPos(1:nPV);
    bestVband = alphaPos(nPV+1:end);
end

function writeLatexAblationTable(S, filename, sysName)
    fid = fopen(filename,'w');
    if fid < 0, error('Could not write %s', filename); end
    fprintf(fid, '%% Auto-generated numeric ablation table for %s\n', sysName);
    fprintf(fid, '\\begin{table}[!t]\n');
    fprintf(fid, '\\caption{Numeric ablation results for the IEEE %s system.}\n', strrep(sysName,'bus','-bus'));
    fprintf(fid, '\\label{tab:ablation_%s}\n', sysName);
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\scriptsize\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Case & Min V & Max V & $J_V$ & Viol. & Red. (\\%%)\\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(S)
        fprintf(fid, '%s & %.4f & %.4f & %.4f & %d & %.1f\\\\\n', ...
            char(S.Case{i}), S.MinV(i), S.MaxV(i), S.TotalVoltageDeviation(i), S.TotalViolationCount(i), S.ViolationReductionPercent(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
