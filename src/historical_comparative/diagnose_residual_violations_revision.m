%% diagnose_residual_violations_revision.m
% Reviewer-revision diagnostic:
% identifies every residual voltage-limit violation in the V4 consistent-Qmax
% "Proposed hybrid" results for the IEEE 33-bus and 69-bus cases.
%
% Run this script from the ALEXIS_PAPER main folder OR from anywhere inside it.
% It searches recursively for the required V4 .mat files and profile files.
%
% Outputs:
%   results_revision_diagnostics/
%       residual_violations_33bus.csv
%       residual_violations_69bus.csv
%       residual_violation_summary.csv
%       residual_violation_report.txt
%
% One "violation" = one bus-time sample outside [0.95, 1.05] p.u.

clear; clc;

Vmin = 0.95;
Vmax = 1.05;
dtHours = 0.25;

rootDir = findProjectRoot(pwd);
outDir = fullfile(rootDir,'results_revision_diagnostics');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('Project root: %s\n', rootDir);

% Preferred files: V4 consistent-Qmax replay, because these are aligned with
% the final component-contribution tables and the corrected inverter-Q model.
f33 = findUniqueRecursive(rootDir, 'ablation_results_33bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');
f69 = findUniqueRecursive(rootDir, 'ablation_results_69bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');

fprintf('33-bus result: %s\n', f33);
fprintf('69-bus result: %s\n', f69);

S33 = load(f33);
S69 = load(f69);

if ~isfield(S33,'result') || ~isfield(S33.result,'Vcase') || ...
        ~isfield(S33.result.Vcase,'ProposedHybrid')
    error('33-bus V4 file does not contain result.Vcase.ProposedHybrid.');
end
if ~isfield(S69,'result') || ~isfield(S69.result,'Vcase') || ...
        ~isfield(S69.result.Vcase,'ProposedHybrid')
    error('69-bus V4 file does not contain result.Vcase.ProposedHybrid.');
end

V33 = S33.result.Vcase.ProposedHybrid;
V69 = S69.result.Vcase.ProposedHybrid;

% Profiles are used only for interpretation columns. The violation locations
% themselves come directly from the saved final voltage matrices.
profile33 = findUniqueRecursive(rootDir, 'myloadspvs.mat', '');
load69file = findUniqueRecursive(rootDir, 'myloads.mat', '');
pv69file   = findUniqueRecursive(rootDir, 'mypvs.mat', '');

P33 = load(profile33);
L69 = load(load69file);
P69 = load(pv69file);

% 33-bus profile preparation: V4 final scaling uses 0.20 MW/site peak.
load33 = P33.myloads(:,1:96);
pv33raw = P33.myPVs(:,1:96);
mx = max(pv33raw(:));
if mx <= 0, error('33-bus PV profile has non-positive maximum.'); end
pv33 = 0.20 * pv33raw / mx;

% 69-bus profiles: final V4 treatment uses uploaded PV values directly in MW.
load69 = L69.myloads(:,1:96);
if isfield(P69,'PV')
    pv69 = P69.PV(:,1:96);
elseif isfield(P69,'mypvs')
    pv69 = P69.mypvs(:,1:96);
else
    error('mypvs.mat does not contain variable PV or mypvs.');
end

pvBuses33 = [21 24 32 17 18];
pvBuses69 = [15 18 27 35 46 60];

T33 = makeViolationTable(V33, load33, pv33, pvBuses33, Vmin, Vmax, dtHours, '33bus');
T69 = makeViolationTable(V69, load69, pv69, pvBuses69, Vmin, Vmax, dtHours, '69bus');

writetable(T33, fullfile(outDir,'residual_violations_33bus.csv'));
writetable(T69, fullfile(outDir,'residual_violations_69bus.csv'));

summary = table();
summary.System = ["33-bus"; "69-bus"];
summary.TotalViolations = [height(T33); height(T69)];
summary.Undervoltage = [sum(T33.Type=="UNDER"); sum(T69.Type=="UNDER")];
summary.Overvoltage = [sum(T33.Type=="OVER"); sum(T69.Type=="OVER")];
summary.MinVoltage = [min(V33(:)); min(V69(:))];
summary.MaxVoltage = [max(V33(:)); max(V69(:))];
summary.UniqueBusesWithViolations = [numel(unique(T33.Bus)); numel(unique(T69.Bus))];
summary.UniqueTimeStepsWithViolations = [numel(unique(T33.TimeStep)); numel(unique(T69.TimeStep))];

writetable(summary, fullfile(outDir,'residual_violation_summary.csv'));

reportFile = fullfile(outDir,'residual_violation_report.txt');
fid = fopen(reportFile,'w');
fprintf(fid,'Residual voltage-violation diagnostic\n');
fprintf(fid,'=====================================\n\n');
fprintf(fid,'Definition: one violation = one bus-time sample with V < %.2f or V > %.2f p.u.\n\n',Vmin,Vmax);

writeSystemReport(fid, T33, V33, '33-bus');
writeSystemReport(fid, T69, V69, '69-bus');

fclose(fid);

disp(summary);

fprintf('\nWorst 33-bus residual samples:\n');
disp(T33(1:min(15,height(T33)),:));

fprintf('\nWorst 69-bus residual samples:\n');
disp(T69(1:min(20,height(T69)),:));

fprintf('\nOutputs saved to:\n%s\n', outDir);
fprintf('Please send residual_violation_report.txt and residual_violation_summary.csv for the manuscript revision.\n');


%% ============================ LOCAL FUNCTIONS ============================

function T = makeViolationTable(V, loadMult, Ppv, pvBuses, Vmin, Vmax, dtHours, systemName)
    [nt, nb] = size(V);
    if nt ~= 96
        warning('%s voltage matrix has %d time steps, expected 96.', systemName, nt);
    end

    [tUnder,bUnder] = find(V < Vmin);
    [tOver,bOver]   = find(V > Vmax);

    t = [tUnder; tOver];
    b = [bUnder; bOver];
    typ = [repmat("UNDER",numel(tUnder),1); repmat("OVER",numel(tOver),1)];

    volt = zeros(numel(t),1);
    deviation = zeros(numel(t),1);
    busLoadMultiplier = nan(numel(t),1);
    totalPV_MW = nan(numel(t),1);
    localPV_MW = zeros(numel(t),1);
    isPVBus = false(numel(t),1);

    for k = 1:numel(t)
        volt(k) = V(t(k),b(k));
        if typ(k) == "UNDER"
            deviation(k) = Vmin - volt(k);
        else
            deviation(k) = volt(k) - Vmax;
        end

        if b(k) <= size(loadMult,1) && t(k) <= size(loadMult,2)
            busLoadMultiplier(k) = loadMult(b(k),t(k));
        end
        if t(k) <= size(Ppv,2)
            totalPV_MW(k) = sum(Ppv(:,t(k)));
        end

        idx = find(pvBuses == b(k),1);
        if ~isempty(idx) && idx <= size(Ppv,1) && t(k) <= size(Ppv,2)
            isPVBus(k) = true;
            localPV_MW(k) = Ppv(idx,t(k));
        end
    end

    Hour = (t-1)*dtHours;
    Clock = strings(numel(t),1);
    for k = 1:numel(t)
        hh = floor(Hour(k));
        mm = round((Hour(k)-hh)*60);
        Clock(k) = sprintf('%02d:%02d',hh,mm);
    end

    T = table(repmat(string(systemName),numel(t),1), t, Hour, Clock, b, typ, ...
              volt, deviation, busLoadMultiplier, totalPV_MW, isPVBus, localPV_MW, ...
        'VariableNames', {'System','TimeStep','Hour','Clock','Bus','Type', ...
        'Voltage_pu','LimitDeviation_pu','BusLoadMultiplier','TotalPV_MW', ...
        'IsPVBus','LocalPV_MW'});

    % Put the most severe residual samples first.
    T = sortrows(T, {'LimitDeviation_pu','TimeStep'}, {'descend','ascend'});
end

function writeSystemReport(fid, T, V, name)
    fprintf(fid,'%s\n',name);
    fprintf(fid,'%s\n',repmat('-',1,length(name)));
    fprintf(fid,'Total residual violations: %d\n',height(T));
    fprintf(fid,'Undervoltage samples: %d\n',sum(T.Type=="UNDER"));
    fprintf(fid,'Overvoltage samples: %d\n',sum(T.Type=="OVER"));
    fprintf(fid,'Minimum voltage: %.6f p.u.\n',min(V(:)));
    fprintf(fid,'Maximum voltage: %.6f p.u.\n',max(V(:)));

    ub = unique(T.Bus);
    ut = unique(T.TimeStep);
    fprintf(fid,'Buses with residual violations: %s\n',mat2str(ub'));
    fprintf(fid,'Time steps with residual violations: %s\n',mat2str(ut'));

    if ~isempty(T)
        fprintf(fid,'Worst sample: bus %d, t=%d (%s), V=%.6f p.u., type=%s\n', ...
            T.Bus(1),T.TimeStep(1),T.Clock(1),T.Voltage_pu(1),T.Type(1));
        fprintf(fid,'Worst 10 samples:\n');
        n = min(10,height(T));
        for k = 1:n
            fprintf(fid,'  %2d) t=%2d (%s), bus=%2d, V=%.6f, loadMult=%.4f, totalPV=%.4f MW, PVbus=%d\n', ...
                k,T.TimeStep(k),T.Clock(k),T.Bus(k),T.Voltage_pu(k), ...
                T.BusLoadMultiplier(k),T.TotalPV_MW(k),T.IsPVBus(k));
        end
    end
    fprintf(fid,'\n');
end

function rootDir = findProjectRoot(startDir)
    % Walk upward until the main scripts or ALEXIS_v4_short_names_bundle appears.
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

function f = findUniqueRecursive(rootDir, fileName, preferredPathToken)
    D = dir(fullfile(rootDir,'**',fileName));
    if isempty(D)
        error('Could not find %s under %s',fileName,rootDir);
    end

    paths = strings(numel(D),1);
    for k = 1:numel(D)
        paths(k) = string(fullfile(D(k).folder,D(k).name));
    end

    if ~isempty(preferredPathToken)
        hit = contains(paths, preferredPathToken, 'IgnoreCase', true);
        if any(hit)
            paths = paths(hit);
        end
    end

    % Prefer files inside ALEXIS_v4_short_names_bundle when duplicates exist.
    hitV4 = contains(paths,'ALEXIS_v4_short_names_bundle','IgnoreCase',true);
    if any(hitV4)
        paths = paths(hitV4);
    end

    f = char(paths(1));
    if numel(paths) > 1
        fprintf('Note: multiple copies of %s found; using:\n  %s\n',fileName,f);
    end
end
