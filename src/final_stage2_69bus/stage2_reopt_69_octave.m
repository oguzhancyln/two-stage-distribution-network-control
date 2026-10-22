function stage2_reopt_69_octave(mode)
% IEEE 69-bus: verify legacy 96-step replay under the 1.00-MW inverter
% rating, then (optionally) re-optimize ONLY the 27 PV-active inverter steps.
%
% USAGE FROM ALEXIS_PAPER PROJECT ROOT:
%   octave --no-gui --quiet --eval "stage2_reopt_69_octave('check')"
%   octave --no-gui --quiet --eval "stage2_reopt_69_octave('optimize')"
%
% IMPORTANT:
%   - No Stage-1 (topology/capacitor) optimization is performed.
%   - The existing source control_schedule_69bus.csv is NEVER modified.
%   - The optimizer is warm-started with the stored inverter parameters,
%     then runs GWO with 40 wolves x 50 iterations for each of 27 steps.
%   - Each step uses a separately recorded deterministic RNG seed.
%   - At each time step the pre-action AC power flow is cached; only
%     post-action AC power flows are used in GWO fitness evaluation.
%   - This is a NEW Stage-2 optimization, not a reconstruction of seed 42.
%
% Copyright: generated for a focused independent reproducibility check.

if nargin < 1 || isempty(mode), mode = 'check'; end
if ~(strcmpi(mode,'check') || strcmpi(mode,'optimize'))
    error('Mode must be ''check'' or ''optimize''.');
end

projectDir = fileparts(mfilename('fullpath'));
bundleDir = fullfile(projectDir, 'ALEXIS_v4_short_names_bundle');
if exist(bundleDir,'dir') ~= 7
    error('Place this file in ALEXIS_PAPER root next to ALEXIS_v4_short_names_bundle.');
end
addpath(bundleDir);
if exist('runpf','file') ~= 2 || exist('loadcase','file') ~= 2
    mpdir = getenv('MATPOWER_DIR');
    if ~isempty(mpdir) && exist(mpdir,'dir') == 7
        addpath(genpath(mpdir));
    end
end
if exist('runpf','file') ~= 2
    error('MATPOWER is not on the Octave path. Use addpath(genpath(''/path/to/matpower'')) or set MATPOWER_DIR.');
end

csvIn = fullfile(bundleDir, 'control_schedule_69bus.csv');
loadFile = fullfile(bundleDir, 'myloads.mat');
pvFile = fullfile(bundleDir, 'mypvs.mat');
assert(exist(csvIn,'file') == 2 && exist(loadFile,'file') == 2 && exist(pvFile,'file') == 2, ...
    'Missing required bundled schedule or load/PV profile.');

cfg = struct();
cfg.pvBuses = [15 18 27 35 46 60];
cfg.capBuses = [7 23 30 39 55];
cfg.capStepMVAr = 0.060;
cfg.PVratedMW = 1.00;
cfg.pfMin = 0.90;
cfg.vMin = 0.95;
cfg.vMax = 1.05;
cfg.penalty = 1e5;
cfg.nWolves = 40;
cfg.nIter = 50;
cfg.seedBase = 690042;
cfg.outputDir = fullfile(projectDir, 'results_stage2_reopt_69_1MW');

mpc0 = loadcase('case69modified');
assert(size(mpc0.bus,1) == 69, 'Wrong number of buses in case69modified.');
assert(size(mpc0.branch,2) >= 11, 'Branch status column is missing.');
Ldata = load(loadFile);
Pdata = load(pvFile);
assert(isfield(Ldata,'myloads') && isfield(Pdata,'PV'), 'Unexpected MAT variable names.');
loadMult = Ldata.myloads;
pvMW = Pdata.PV;
assert(isequal(size(loadMult), [69 96]) && isequal(size(pvMW), [6 96]), ...
    'Profiles do not have expected dimensions (69x96 loads, 6x96 PV).');
assert(all(isfinite(loadMult(:))) && all(isfinite(pvMW(:))), 'Nonfinite profile values.');
assert(all(pvMW(:) >= -1e-12) && max(pvMW(:)) <= cfg.PVratedMW + 1e-9, ...
    'PV output is incompatible with 1.00-MW rated inverter.');

schedule = s2_read_schedule(csvIn);
activeTimes = [25:39 61:72];
assert(numel(activeTimes) == 27, 'Unexpected active-window length.');
for t = 1:96
    assert(schedule(t).t == t, 'Schedule time step order error.');
    assert(numel(schedule(t).open) == size(mpc0.branch,1) - 68, ...
        'Incorrect number of open lines in schedule.');
    assert(numel(unique(schedule(t).open)) == numel(schedule(t).open), ...
        'Duplicate open lines in schedule.');
    assert(numel(schedule(t).caps) == 5, 'Incorrect number of capacitor entries.');
    if any(activeTimes == t)
        assert(strcmp(schedule(t).mode, 'DROOP') && ...
            numel(schedule(t).vsp) == 6 && numel(schedule(t).vband) == 6, ...
            'Expected DROOP row missing or malformed.');
    else
        assert(strcmp(schedule(t).mode, 'DNR+CAP'), ...
            'Unexpected control mode outside inverter windows.');
    end
end

mpopt = mpoption('verbose',0,'out.all',0);
legacyV = zeros(96,69);
legacyScore = zeros(96,1);
preV = zeros(96,69);
legacyAccepted = false(96,1);
legacyQ = zeros(96,6);
prepared = cell(96,1);

fprintf('Preflight: replaying all 96 time steps with 1.00 MW/site and pf_min = 0.90...\n');
for t = 1:96
    prepared{t} = s2_prepare(mpc0, loadMult, pvMW, schedule(t), t, cfg);
    [preV(t,:), ok] = s2_pf(prepared{t}, mpopt);
    if ~ok, error('Pre-action PF failed at t=%d. No optimization executed.', t); end
    if any(activeTimes == t)
        [postV, postScore, ok, q] = s2_post(prepared{t}, preV(t,:), ...
            schedule(t).vsp, schedule(t).vband, pvMW(:,t), cfg, mpopt);
        preScore = s2_obj(preV(t,:),cfg);
        if ok && postScore <= preScore
            legacyV(t,:) = postV;
            legacyScore(t) = postScore;
            legacyQ(t,:) = q(:)';
            legacyAccepted(t) = true;
        else
            legacyV(t,:) = preV(t,:);
            legacyScore(t) = preScore;
        end
    else
        legacyV(t,:) = preV(t,:);
        legacyScore(t) = s2_obj(preV(t,:),cfg);
    end
end
legacy = s2_metrics(legacyV,cfg);
fprintf('Legacy schedule, 1.00-MW replay: JV=%.12f, violations=%d, Vmin=%.9f, Vmax=%.9f\n', ...
    legacy.JV, legacy.viol, legacy.minV, legacy.maxV);
fprintf('Old droop commands accepted by replayer: %d/27\n',nnz(legacyAccepted));

% Published V4 replay (saved in ablation_summary_69bus.csv):
expectedJV = 1.58073212457191;
expectedViol = 51;
if abs(legacy.JV - expectedJV) > 1e-4 || legacy.viol ~= expectedViol
    error(['Preflight does not reproduce the published V4 result ' ...
      '(expected JV=%.12f, violations=%d). STOP. Send terminal log; ' ...
      'no optimizer was launched.'], expectedJV, expectedViol);
end
fprintf('PASS: exact result family matches the published V4 69-bus replay.\n');
if strcmpi(mode,'check')
    fprintf('CHECK mode complete: no GWO optimization started; no source files modified.\n');
    return;
end

if exist(cfg.outputDir,'dir') ~= 7, mkdir(cfg.outputDir); end
finalCsv = fullfile(cfg.outputDir,'control_schedule_69bus_STAGE2_REOPT.csv');
if exist(finalCsv,'file') == 2
    error('Final output already exists; refusing to overwrite: %s', finalCsv);
end

newSchedule = schedule;
newV = preV;
newQ = zeros(96,6);
accepted = false(96,1);
fitnessCalls = 0;
allCurves = nan(numel(activeTimes),cfg.nIter);

fprintf('\nRe-optimizing only 27 inverter time steps (40 wolves x 50 iterations).\n');
fprintf('Saved DNR and CAP states will not be changed.\n');
for k = 1:numel(activeTimes)
    t = activeTimes(k);
    checkpoint = fullfile(cfg.outputDir,sprintf('step_%02d.mat',t));
    if exist(checkpoint,'file') == 2
        loaded = load(checkpoint);
        assert(isfield(loaded,'vspNew') && isfield(loaded,'vbandNew') && ...
               isfield(loaded,'curve') && isfield(loaded,'localCalls'), ...
               'Invalid checkpoint format at t=%d.',t);
        vspNew = loaded.vspNew;
        vbandNew = loaded.vbandNew;
        curve = loaded.curve;
        localCalls = loaded.localCalls;
        fprintf('t=%02d (%02d/27): loaded previously saved checkpoint.\n', t,k);
    else
        rng(cfg.seedBase + t,'twister');
        oldX = [schedule(t).vsp schedule(t).vband];
        [vspNew,vbandNew,curve,localCalls] = s2_gwo(prepared{t},preV(t,:), ...
            pvMW(:,t),oldX,cfg,mpopt,t,k);
        save(checkpoint,'vspNew','vbandNew','curve','localCalls','t');
    end
    fitnessCalls = fitnessCalls + localCalls;
    allCurves(k,:) = curve;
    [postV,postScore,ok,q] = s2_post(prepared{t},preV(t,:), ...
        vspNew,vbandNew,pvMW(:,t),cfg,mpopt);
    preScore = s2_obj(preV(t,:),cfg);
    if ok && postScore <= preScore
        newV(t,:) = postV;
        newQ(t,:) = q(:)';
        newSchedule(t).mode = 'DROOP';
        newSchedule(t).vsp = vspNew;
        newSchedule(t).vband = vbandNew;
        newSchedule(t).score = postScore;
        accepted(t) = true;
    else
        newV(t,:) = preV(t,:);
        newSchedule(t).mode = 'DNR+CAP-DROOP-NOT-USED';
        newSchedule(t).vsp = [];
        newSchedule(t).vband = [];
        newSchedule(t).score = preScore;
    end
end
newMetrics = s2_metrics(newV,cfg);
fprintf('\nOLD 1MW REPLAY: JV=%.12f, violations=%d, Vmin=%.9f\n', ...
    legacy.JV,legacy.viol,legacy.minV);
fprintf('NEW STAGE2 REOPT: JV=%.12f, violations=%d, Vmin=%.9f\n', ...
    newMetrics.JV,newMetrics.viol,newMetrics.minV);
fprintf('Inverter steps accepted: %d/27; total GWO fitness calls: %d.\n', ...
    nnz(accepted),fitnessCalls);

s2_write_schedule(finalCsv,newSchedule);
s2_write_summary(fullfile(cfg.outputDir,'comparison_69bus.csv'),legacy,newMetrics);
save(fullfile(cfg.outputDir,'stage2_reopt_69_results.mat'), ...
    'newSchedule','newV','preV','legacyV','legacy','newMetrics', ...
    'newQ','accepted','allCurves','cfg','activeTimes','fitnessCalls');
fprintf('Wrote NEW schedule and comparison to %s\n',cfg.outputDir);
fprintf('Original schedule and all existing published results were left untouched.\n');
end

function schedule = s2_read_schedule(filename)
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot read schedule CSV: %s',filename); end
    header = fgetl(fid);
    if ~ischar(header) || isempty(strfind(header,'TimeStep'))
        fclose(fid); error('Unexpected schedule CSV header.');
    end
    schedule = repmat(struct('t',0,'mode','','open',[],'caps',[], ...
        'vsp',[],'vband',[],'score',NaN),96,1);
    for i=1:96
        line = fgetl(fid);
        if ~ischar(line), fclose(fid); error('Schedule CSV ended before row 96.'); end
        f = strsplit(line,',');
        if numel(f) ~= 7, fclose(fid); error('CSV row %d has %d fields, expected 7.',i,numel(f)); end
        schedule(i).t = str2double(f{1});
        schedule(i).mode = strtrim(f{2});
        schedule(i).open = s2_vector(f{3});
        schedule(i).caps = s2_vector(f{4});
        schedule(i).vsp = s2_vector(f{5});
        schedule(i).vband = s2_vector(f{6});
        schedule(i).score = str2double(f{7});
    end
    fclose(fid);
end

function v = s2_vector(s)
    s = strtrim(strrep(strrep(s,'[',''),']',''));
    if isempty(s), v = []; return; end
    v = sscanf(s,'%f').';
end

function m = s2_prepare(mpc,loads,pv,sched,t,cfg)
    m = mpc;
    m.branch(:,11) = 1;
    m.branch(sched.open,11) = 0;
    m.bus(:,3) = mpc.bus(:,3).*loads(:,t);
    m.bus(:,4) = mpc.bus(:,4).*loads(:,t);
    for j=1:6
        b = cfg.pvBuses(j);
        m.bus(b,3) = m.bus(b,3)-pv(j,t);
    end
    for j=1:5
        b = cfg.capBuses(j);
        m.bus(b,4) = m.bus(b,4)-sched.caps(j)*cfg.capStepMVAr;
    end
end

function [v,ok] = s2_pf(m,opts)
    try
        r = runpf(m,opts);
        ok = isfield(r,'success') && (r.success == 1);
        if ok
            v = abs(r.bus(:,8)).';
            ok = all(isfinite(v));
        else
            v = nan(1,size(m.bus,1));
        end
    catch
        v = nan(1,size(m.bus,1));
        ok = false;
    end
end

function [v,J,ok,q] = s2_post(prepared,vpre,vsp,vband,pvMW,cfg,opts)
    ratedS = cfg.PVratedMW/cfg.pfMin;
    qMax = min(sqrt(max(ratedS^2-pvMW(:).^2,0)), ...
        pvMW(:)*tan(acos(cfg.pfMin)));
    q = zeros(6,1);
    m = prepared;
    for j=1:6
        if qMax(j) <= 0
            q(j) = 0;
        else
            dv = vpre(cfg.pvBuses(j))-vsp(j);
            q(j) = max(-qMax(j),min(qMax(j),-qMax(j)*dv/vband(j)));
        end
        m.bus(cfg.pvBuses(j),4) = m.bus(cfg.pvBuses(j),4)-q(j);
    end
    [v,ok] = s2_pf(m,opts);
    if ok, J = s2_obj(v,cfg); else, J = 1e12; end
end

function J = s2_obj(v,cfg)
    if any(~isfinite(v)), J=1e12; return; end
    J = sum((v-1).^2) + cfg.penalty * ...
        sum(max(cfg.vMin-v,0).^2 + max(v-cfg.vMax,0).^2);
end

function metric = s2_metrics(v,cfg)
    metric.JV = sum((v(:)-1).^2);
    metric.viol = nnz(v(:)<cfg.vMin | v(:)>cfg.vMax);
    metric.minV = min(v(:));
    metric.maxV = max(v(:));
    metric.meanV = mean(v(:));
    metric.under = nnz(v(:)<cfg.vMin);
    metric.over = nnz(v(:)>cfg.vMax);
end

function [vsp,vband,curve,calls] = s2_gwo(prepared,vpre,pvMW,oldX,cfg,opts,t,index)
    dim = 12;
    lb = [0.95*ones(1,6) 0.01*ones(1,6)];
    ub = [1.05*ones(1,6) 0.05*ones(1,6)];
    X = repmat(lb,cfg.nWolves,1) + rand(cfg.nWolves,dim).*repmat(ub-lb,cfg.nWolves,1);
    X(1,:) = min(ub,max(lb,oldX));  % documented warm start
    alphaX = X(1,:); betaX = alphaX; deltaX = alphaX;
    alphaScore = inf; betaScore=inf; deltaScore=inf;
    calls = 0; curve = nan(1,cfg.nIter);
    tstart = tic;
    for iter = 1:cfg.nIter
        for i=1:cfg.nWolves
            X(i,:) = s2_reflect(X(i,:),lb,ub);
            [~,f,ok,~] = s2_post(prepared,vpre,X(i,1:6),X(i,7:12),pvMW,cfg,opts);
            if ~ok, f=1e12; end
            calls = calls+1;
            if f < alphaScore
                deltaScore=betaScore; deltaX=betaX;
                betaScore=alphaScore; betaX=alphaX;
                alphaScore=f; alphaX=X(i,:);
            elseif f < betaScore
                deltaScore=betaScore; deltaX=betaX;
                betaScore=f; betaX=X(i,:);
            elseif f < deltaScore
                deltaScore=f; deltaX=X(i,:);
            end
        end
        a = 2 - iter*(2/cfg.nIter);
        for i=1:cfg.nWolves
            for j=1:dim
                r1=rand();r2=rand();A1=2*a*r1-a;C1=2*r2;
                X1=alphaX(j)-A1*abs(C1*alphaX(j)-X(i,j));
                r1=rand();r2=rand();A2=2*a*r1-a;C2=2*r2;
                X2=betaX(j)-A2*abs(C2*betaX(j)-X(i,j));
                r1=rand();r2=rand();A3=2*a*r1-a;C3=2*r2;
                X3=deltaX(j)-A3*abs(C3*deltaX(j)-X(i,j));
                X(i,j)=(X1+X2+X3)/3;
            end
            X(i,:) = s2_reflect(X(i,:),lb,ub);
        end
        curve(iter)=alphaScore;
        if mod(iter,10)==0 || iter==1 || iter==cfg.nIter
            fprintf('  t=%02d (%d/27), iter=%d/%d, Jpen=%.8g, elapsed=%.1fs\n', ...
                t,index,iter,cfg.nIter,alphaScore,toc(tstart));
            fflush(stdout);
        end
    end
    vsp = alphaX(1:6);
    vband = alphaX(7:12);
end

function x = s2_reflect(x,lb,ub)
    for j=1:numel(x)
        d = ub(j)-lb(j);
        x(j)=lb(j)+abs(mod(x(j)-ub(j),2*d)-d);
        x(j)=min(ub(j),max(lb(j),x(j)));
    end
end

function s2_write_schedule(filename,rows)
    fid=fopen(filename,'w');
    if fid<0, error('Cannot write schedule CSV.'); end
    fprintf(fid,'TimeStep,Mode,OpenLines,CapSteps,Vsp,Vband,Score\n');
    for t=1:96
        fprintf(fid,'%d,%s,%s,%s,%s,%s,%.15g\n', t,rows(t).mode, ...
            mat2str(rows(t).open),mat2str(rows(t).caps), ...
            mat2str(rows(t).vsp,12),mat2str(rows(t).vband,12),rows(t).score);
    end
    fclose(fid);
end

function s2_write_summary(filename,old,new)
    fid=fopen(filename,'w');
    if fid<0, error('Cannot write comparison.'); end
    fprintf(fid,'Method,JV,Violations,MinV,MaxV,MeanV,UnderVoltage,OverVoltage\n');
    fprintf(fid,'Legacy schedule re-evaluated with 1MW,%.15g,%d,%.15g,%.15g,%.15g,%d,%d\n', ...
        old.JV,old.viol,old.minV,old.maxV,old.meanV,old.under,old.over);
    fprintf(fid,'Retained Stage1 plus reoptimized Stage2 at 1MW,%.15g,%d,%.15g,%.15g,%.15g,%d,%d\n', ...
        new.JV,new.viol,new.minV,new.maxV,new.meanV,new.under,new.over);
    fclose(fid);
end
