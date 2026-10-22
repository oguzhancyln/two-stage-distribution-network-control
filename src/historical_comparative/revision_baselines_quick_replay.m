%% revision_baselines_quick_replay.m
% Quick reviewer-revision baselines using the FINAL V4 consistent-Qmax schedules.
%
% Purpose:
%   Add low-cost ablations/policy baselines without re-running the expensive
%   topology optimization. These are "component-removal / policy-replay"
%   diagnostics around the already selected DNR+CAP schedule.
%
% Cases evaluated for BOTH IEEE 33-bus and 69-bus systems:
%   1) No control
%   2) DNR only           : selected topology schedule, capacitors disabled
%   3) CAP only           : default radial topology, selected cap schedule
%   4) DNR+CAP only       : selected topology + capacitor schedule
%   5) Fixed VVC          : DNR+CAP + fixed Vsp=1.00, Vband=0.03 in PV-active steps
%   6) Full no gate       : DNR+CAP + stored optimized droop, always applied
%   7) Proposed hybrid    : DNR+CAP + stored optimized droop + acceptance gate
%   8) Fixed-Q hybrid     : same as proposed, but with a rated-Q limit instead
%                           of the adopted P-dependent PF limit
%
% IMPORTANT:
%   DNR-only and CAP-only here are component-removal ablations of the final
%   selected DNR+CAP schedule. They are NOT separately re-optimized stand-alone
%   controllers. State this explicitly in the paper if these results are used.
%
% Q-limit definitions:
%   Dynamic (paper):
%     Qmax = min( sqrt(S^2-P^2), P*tan(acos(PFmin)) )
%
%   Fixed-rated comparison:
%     Qmax = min( sqrt(S^2-P^2), Prated*tan(acos(PFmin)) )
%   This preserves the inverter apparent-power circle but removes the
%   P-dependent PF-policy shrinkage at low active power.
%
% Outputs:
%   results_revision_baselines/
%       revision_baselines_33bus.csv
%       revision_baselines_69bus.csv
%       revision_baselines_combined.csv
%       revision_baselines_report.txt
%       revision_baselines_33bus.mat
%       revision_baselines_69bus.mat
%
% Run with MATPOWER and case69modified on the MATLAB path.

clear; clc;
rng(42);

cfg.T = 96;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambdaUnder = 1e5;
cfg.lambdaOver = 1e5;
cfg.PFmin = 0.90;
cfg.fixedVsp = 1.00;
cfg.fixedVband = 0.03;

rootDir = findProjectRoot(pwd);
outDir = fullfile(rootDir,'results_revision_baselines');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('Project root: %s\n',rootDir);
fprintf('Output dir  : %s\n\n',outDir);

% Use the V4 consistent-Qmax result files because they contain the corrected
% schedules and system definitions used for the final component analysis.
f33 = findUniqueRecursive(rootDir,'ablation_results_33bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');
f69 = findUniqueRecursive(rootDir,'ablation_results_69bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');

A33 = load(f33);
A69 = load(f69);

if ~isfield(A33,'result') || ~isfield(A33,'sys')
    error('33-bus V4 file must contain result and sys.');
end
if ~isfield(A69,'result') || ~isfield(A69,'sys')
    error('69-bus V4 file must contain result and sys.');
end

sys33 = A33.sys;
sys69 = A69.sys;

% Force the final paper-consistent settings in case an older saved sys differs.
sys33.caseName = 'case33bw';
sys33.pvBuses = [21 24 32 17 18];
sys33.capBuses = [18 25 33];
sys33.capStepMVAr = 60/1000;
sys33.pvRatedMW = 0.20;
sys33.name = '33bus';

sys69.caseName = 'case69modified';
sys69.pvBuses = [15 18 27 35 46 60];
sys69.capBuses = [7 23 30 39 55];
sys69.capStepMVAr = 60/1000;
sys69.pvRatedMW = 1.00;
sys69.name = '69bus';

% Profiles
profile33 = findUniqueRecursive(rootDir,'myloadspvs.mat','');
load69file = findUniqueRecursive(rootDir,'myloads.mat','');
pv69file   = findUniqueRecursive(rootDir,'mypvs.mat','');

P33 = load(profile33);
L69 = load(load69file);
P69 = load(pv69file);

load33 = P33.myloads(:,1:cfg.T);
pv33raw = P33.myPVs(:,1:cfg.T);
mx = max(pv33raw(:));
if mx <= 0, error('33-bus PV profile maximum must be positive.'); end
pv33 = sys33.pvRatedMW * pv33raw / mx;

load69 = L69.myloads(:,1:cfg.T);
if isfield(P69,'PV')
    pv69 = P69.PV(:,1:cfg.T);
elseif isfield(P69,'mypvs')
    pv69 = P69.mypvs(:,1:cfg.T);
else
    error('mypvs.mat must contain PV or mypvs.');
end

fprintf('33-bus PV peak used: %.4f MW/site\n',max(pv33(:)));
fprintf('69-bus PV peak used: %.4f MW/site\n\n',max(pv69(:)));

res33 = runOneSystem(sys33, A33.result.schedule, load33, pv33, cfg);
res69 = runOneSystem(sys69, A69.result.schedule, load69, pv69, cfg);

writetable(res33.summary, fullfile(outDir,'revision_baselines_33bus.csv'));
writetable(res69.summary, fullfile(outDir,'revision_baselines_69bus.csv'));

combined = [res33.summary; res69.summary];
writetable(combined, fullfile(outDir,'revision_baselines_combined.csv'));

save(fullfile(outDir,'revision_baselines_33bus.mat'),'res33','sys33','cfg');
save(fullfile(outDir,'revision_baselines_69bus.mat'),'res69','sys69','cfg');

% Human-readable report
fid = fopen(fullfile(outDir,'revision_baselines_report.txt'),'w');
fprintf(fid,'Reviewer-revision quick baseline report\n');
fprintf(fid,'=======================================\n\n');
fprintf(fid,'NOTE: DNR-only and CAP-only are component-removal replays of the selected\n');
fprintf(fid,'DNR+CAP schedule; they are not separately re-optimized stand-alone controls.\n\n');
writeReport(fid,res33.summary,'33-bus');
writeReport(fid,res69.summary,'69-bus');
fclose(fid);

disp(res33.summary);
disp(res69.summary);

% Sanity check against known final V4 proposed-hybrid headline values.
checkHeadline(res33.summary,'33bus',1.5771,16);
checkHeadline(res69.summary,'69bus',1.5807,51);

fprintf('\nDone. Please send:\n');
fprintf('  %s\n',fullfile(outDir,'revision_baselines_report.txt'));
fprintf('  %s\n',fullfile(outDir,'revision_baselines_combined.csv'));


%% ========================================================================
function out = runOneSystem(sys, sched, loadMult, Ppv, cfg)

    mpc0 = loadcase(sys.caseName);
    nb = size(mpc0.bus,1);
    nbr = size(mpc0.branch,1);
    nOpen = nbr - (nb-1);
    baseOpen = determineDefaultOpenLines(mpc0,nOpen);

    cases = { ...
        'No control', ...
        'DNR only (schedule ablation)', ...
        'CAP only (schedule ablation)', ...
        'DNR+CAP only', ...
        'Fixed VVC on DNR+CAP', ...
        'Full optimized droop - no gate', ...
        'Proposed hybrid', ...
        'Fixed-Q capability hybrid'};

    Vcases = struct();
    accRej = zeros(numel(cases),2);

    for c = 1:numel(cases)
        name = cases{c};
        key = matlab.lang.makeValidName(name);
        fprintf('\n[%s] %s\n',sys.name,name);

        Vmat = nan(cfg.T,nb);
        nAccepted = 0;
        nRejected = 0;

        for t = 1:cfg.T
            row = sched(t);

            % Defaults
            openLines = baseOpen;
            capSteps = zeros(1,numel(sys.capBuses));
            useDroop = false;
            Vsp = [];
            Vband = [];
            gate = true;
            qMode = 'dynamic_pf';

            switch name
                case 'No control'
                    % defaults

                case 'DNR only (schedule ablation)'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end

                case 'CAP only (schedule ablation)'
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end

                case 'DNR+CAP only'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end

                case 'Fixed VVC on DNR+CAP'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end
                    if row.HasDroop
                        useDroop = true;
                        Vsp = cfg.fixedVsp*ones(1,numel(sys.pvBuses));
                        Vband = cfg.fixedVband*ones(1,numel(sys.pvBuses));
                    end

                case 'Full optimized droop - no gate'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end
                    if row.HasDroop
                        useDroop = true;
                        Vsp = row.Vsp;
                        Vband = row.Vband;
                        gate = false;
                    end

                case 'Proposed hybrid'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end
                    if row.HasDroop
                        useDroop = true;
                        Vsp = row.Vsp;
                        Vband = row.Vband;
                    end

                case 'Fixed-Q capability hybrid'
                    if ~isempty(row.OpenLines), openLines = row.OpenLines; end
                    if ~isempty(row.CapSteps), capSteps = row.CapSteps; end
                    if row.HasDroop
                        useDroop = true;
                        Vsp = row.Vsp;
                        Vband = row.Vband;
                        qMode = 'fixed_rated';
                    end
            end

            [V,~,accepted,rejected] = simulateTimeStep(mpc0,openLines,capSteps, ...
                useDroop,Vsp,Vband,t,loadMult,Ppv,sys,cfg,gate,qMode);

            Vmat(t,:) = V;
            nAccepted = nAccepted + accepted;
            nRejected = nRejected + rejected;
        end

        Vcases.(key) = Vmat;
        accRej(c,:) = [nAccepted nRejected];
    end

    S = table();
    S.System = repmat(string(sys.name),numel(cases),1);
    S.Case = string(cases(:));

    for c = 1:numel(cases)
        key = matlab.lang.makeValidName(cases{c});
        met = computeMetrics(Vcases.(key),cfg);
        S.MinV(c,1) = met.MinV;
        S.MaxV(c,1) = met.MaxV;
        S.MeanV(c,1) = met.MeanV;
        S.JV(c,1) = met.TotalVoltageDeviation;
        S.UnderVoltageCount(c,1) = met.UnderVoltageCount;
        S.OverVoltageCount(c,1) = met.OverVoltageCount;
        S.TotalViolationCount(c,1) = met.TotalViolationCount;
        S.InverterAccepted(c,1) = accRej(c,1);
        S.InverterRejected(c,1) = accRej(c,2);
    end

    baseJ = S.JV(1);
    baseViol = S.TotalViolationCount(1);
    S.JVReductionPercent = 100*(baseJ-S.JV)./max(baseJ,eps);
    S.ViolationReductionPercent = 100*(baseViol-S.TotalViolationCount)./max(baseViol,eps);

    out.summary = S;
    out.Vcases = Vcases;
    out.baseOpen = baseOpen;
end


function [Vfinal,score,accepted,rejected] = simulateTimeStep( ...
    mpc0,openLines,capSteps,useDroop,Vsp,Vband,t,loadMult,Ppv,sys,cfg,gate,qMode)

    accepted = 0;
    rejected = 0;

    m0 = prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys);
    mpopt = mpoption('verbose',0,'out.all',0);

    if ~useDroop
        r = runpf(m0,mpopt);
        if ~r.success
            Vfinal = nan(1,size(mpc0.bus,1));
            score = inf;
            return;
        end
        Vfinal = abs(r.bus(:,8)).';
        score = objectiveVoltage(Vfinal,cfg);
        return;
    end

    r0 = runpf(m0,mpopt);
    if ~r0.success
        Vfinal = nan(1,size(mpc0.bus,1));
        score = inf;
        return;
    end

    V0 = abs(r0.bus(:,8)).';
    score0 = objectiveVoltage(V0,cfg);

    QmaxStep = computeQmax(Ppv(:,t),mpc0.baseMVA,cfg,sys,qMode);

    m1 = m0;
    for k = 1:numel(sys.pvBuses)
        b = sys.pvBuses(k);
        Qpu = pvControlSetPointLocal(V0(b),Vsp(k),Vband(k),QmaxStep(k));
        m1.bus(b,4) = m1.bus(b,4) - Qpu*mpc0.baseMVA;
    end

    r1 = runpf(m1,mpopt);
    if ~r1.success
        Vfinal = V0;
        score = score0;
        rejected = 1;
        return;
    end

    V1 = abs(r1.bus(:,8)).';
    score1 = objectiveVoltage(V1,cfg);

    if gate && score1 > score0
        Vfinal = V0;
        score = score0;
        rejected = 1;
    else
        Vfinal = V1;
        score = score1;
        accepted = 1;
    end
end


function m = prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys)
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

    if ~isempty(capSteps)
        for k = 1:min(numel(sys.capBuses),numel(capSteps))
            m.bus(sys.capBuses(k),4) = m.bus(sys.capBuses(k),4) ...
                - capSteps(k)*sys.capStepMVAr;
        end
    end
end


function QmaxStep = computeQmax(PactMW,baseMVA,cfg,sys,qMode)
    SmaxMW = sys.pvRatedMW ./ cfg.PFmin;
    qCircleMVAr = sqrt(max(SmaxMW.^2 - PactMW.^2,0));

    switch qMode
        case 'dynamic_pf'
            qPolicyMVAr = PactMW .* tan(acos(cfg.PFmin));

        case 'fixed_rated'
            qRatedMVAr = sys.pvRatedMW * tan(acos(cfg.PFmin));
            qPolicyMVAr = qRatedMVAr * ones(size(PactMW));

        otherwise
            error('Unknown qMode: %s',qMode);
    end

    qMaxMVAr = min(qCircleMVAr,qPolicyMVAr);
    QmaxStep = qMaxMVAr ./ baseMVA;
end


function Qcmd = pvControlSetPointLocal(Vmeas,Vsp,Vband,Qmax)
    if Qmax <= 0 || Vband <= 0
        Qcmd = 0;
        return;
    end
    K = Qmax ./ Vband;
    dV = Vmeas - Vsp;
    Qcmd = -K .* dV;
    if dV <= -Vband
        Qcmd = Qmax;
    elseif dV >= Vband
        Qcmd = -Qmax;
    end
    Qcmd = min(Qmax,max(-Qmax,Qcmd));
end


function J = objectiveVoltage(V,cfg)
    if any(isnan(V)) || any(isinf(V))
        J = inf;
        return;
    end
    under = max(cfg.Vmin-V,0);
    over = max(V-cfg.Vmax,0);
    J = sum((V-1).^2) ...
        + cfg.lambdaUnder*sum(under.^2) ...
        + cfg.lambdaOver*sum(over.^2);
end


function met = computeMetrics(Vmat,cfg)
    V = Vmat(:);
    V = V(~isnan(V));
    met.MinV = min(V);
    met.MaxV = max(V);
    met.MeanV = mean(V);
    met.TotalVoltageDeviation = sum((V-1).^2);
    met.UnderVoltageCount = sum(V<cfg.Vmin);
    met.OverVoltageCount = sum(V>cfg.Vmax);
    met.TotalViolationCount = met.UnderVoltageCount+met.OverVoltageCount;
end


function openLines = determineDefaultOpenLines(mpc,nOpen)
    existingOpen = find(mpc.branch(:,11)==0).';
    if numel(existingOpen)==nOpen
        openLines = existingOpen;
    else
        openLines = (size(mpc.branch,1)-nOpen+1):size(mpc.branch,1);
    end
end


function writeReport(fid,S,name)
    fprintf(fid,'%s\n',name);
    fprintf(fid,'%s\n',repmat('-',1,length(name)));
    for i = 1:height(S)
        fprintf(fid,'%-34s MinV=%8.5f MaxV=%8.5f JV=%9.4f Viol=%4d Acc=%3d Rej=%3d\n', ...
            char(S.Case(i)),S.MinV(i),S.MaxV(i),S.JV(i), ...
            S.TotalViolationCount(i),S.InverterAccepted(i),S.InverterRejected(i));
    end
    fprintf(fid,'\n');
end


function checkHeadline(S,systemName,targetJ,targetViol)
    idx = find(S.Case=="Proposed hybrid",1);
    if isempty(idx), return; end
    dJ = abs(S.JV(idx)-targetJ);
    dV = abs(S.TotalViolationCount(idx)-targetViol);
    if dJ > 5e-3 || dV > 0
        warning('%s replay does not exactly match headline: JV %.4f vs %.4f; violations %d vs %d.', ...
            systemName,S.JV(idx),targetJ,S.TotalViolationCount(idx),targetViol);
    else
        fprintf('%s proposed-hybrid replay matches headline closely.\n',systemName);
    end
end


function rootDir = findProjectRoot(startDir)
    d = startDir;
    while true
        if exist(fullfile(d,'ALEXIS_v4_short_names_bundle'),'dir') || ...
           exist(fullfile(d,'hybrid_dnr_cap_droop_33bus_final_v2_safe.m'),'file')
            rootDir = d;
            return;
        end
        parent = fileparts(d);
        if strcmp(parent,d)
            rootDir = startDir;
            return;
        end
        d = parent;
    end
end


function f = findUniqueRecursive(rootDir,fileName,preferredToken)
    D = dir(fullfile(rootDir,'**',fileName));
    if isempty(D)
        error('Could not find %s under %s',fileName,rootDir);
    end

    paths = strings(numel(D),1);
    for k = 1:numel(D)
        paths(k) = string(fullfile(D(k).folder,D(k).name));
    end

    if ~isempty(preferredToken)
        hit = contains(paths,preferredToken,'IgnoreCase',true);
        if any(hit), paths = paths(hit); end
    end

    hitV4 = contains(paths,'ALEXIS_v4_short_names_bundle','IgnoreCase',true);
    if any(hitV4), paths = paths(hitV4); end

    f = char(paths(1));
    if numel(paths)>1
        fprintf('Multiple copies of %s found; using:\n  %s\n',fileName,f);
    end
end
