%% revision_operational_feasibility_diagnostics.m
% Reviewer-revision diagnostic for the FINAL two-stage schedules.
%
% This script does NOT add new constraints to the optimizer. It checks what
% can be verified ex post from the saved final schedules and MATPOWER case data:
%   - AC power-flow convergence
%   - connected radial endpoint topology
%   - hard voltage-limit satisfaction
%   - thermal loading IF positive RATE_A values exist in the case file
%   - substation reverse active-power flow
%
% It also reports limitations that cannot be verified from the current model:
%   - no explicit field switchability set
%   - no intermediate switching-sequence feasibility/protection model
%
% Outputs:
%   results_revision_operational_feasibility/
%       operational_feasibility_summary.csv
%       operational_feasibility_report.txt
%
% MATPOWER and case69modified must be on the MATLAB path.

clear; clc;
rng(42);

cfg.T = 96;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.PFmin = 0.90;

rootDir = findProjectRoot(pwd);
outDir = fullfile(rootDir,'results_revision_operational_feasibility');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('Project root: %s\n',rootDir);

f33 = findUniqueRecursive(rootDir,'ablation_results_33bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');
f69 = findUniqueRecursive(rootDir,'ablation_results_69bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');

A33 = load(f33);
A69 = load(f69);

sys33 = A33.sys;
sys69 = A69.sys;

% Final paper-consistent settings
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

R33 = inspectSystem(sys33,A33.result.schedule,load33,pv33,cfg);
R69 = inspectSystem(sys69,A69.result.schedule,load69,pv69,cfg);

S = [R33.summary; R69.summary];
writetable(S,fullfile(outDir,'operational_feasibility_summary.csv'));

fid = fopen(fullfile(outDir,'operational_feasibility_report.txt'),'w');
fprintf(fid,'Operational-feasibility diagnostic for selected two-stage schedules\n');
fprintf(fid,'================================================================\n\n');
fprintf(fid,'Scope: ex-post diagnostic only. The optimization itself does not enforce\n');
fprintf(fid,'hard voltage limits, thermal limits, reverse-flow limits, or switching-sequence constraints.\n\n');
writeReport(fid,R33);
writeReport(fid,R69);

fprintf(fid,'Model limitations common to both feeders\n');
fprintf(fid,'-----------------------------------------\n');
fprintf(fid,'1) Voltage limits are soft objective penalties, not hard constraints.\n');
fprintf(fid,'2) Thermal feasibility is testable only if positive MATPOWER RATE_A values exist.\n');
fprintf(fid,'3) Reverse flow is not constrained by the optimization; this report only checks slack active power ex post.\n');
fprintf(fid,'4) The implementation does not define a restricted set of physically switchable branches;\n');
fprintf(fid,'   decoded open-branch decisions are treated as controllable branch statuses in the benchmark model.\n');
fprintf(fid,'5) Only endpoint radial topologies are validated. Intermediate switching sequences,\n');
fprintf(fid,'   transient make-before-break/break-before-make states, protection coordination,\n');
fprintf(fid,'   switching times, and communication delays are outside the model.\n');
fclose(fid);

disp(S);

fprintf('\nDone. Please send:\n');
fprintf('  %s\n',fullfile(outDir,'operational_feasibility_report.txt'));
fprintf('  %s\n',fullfile(outDir,'operational_feasibility_summary.csv'));


%% ========================================================================
function R = inspectSystem(sys,sched,loadMult,Ppv,cfg)

    mpc0 = loadcase(sys.caseName);
    nb = size(mpc0.bus,1);
    nbr = size(mpc0.branch,1);
    nOpen = nbr-(nb-1);

    % MATPOWER indices
    F_BUS=1; T_BUS=2; RATE_A=6; BR_STATUS=11;
    PG=2; GEN_BUS=1;
    PF=14; QF=15; PT=16; QT=17;

    mpopt = mpoption('verbose',0,'out.all',0);

    pfConverged = false(cfg.T,1);
    radialOK = false(cfg.T,1);
    minV = nan(cfg.T,1);
    maxV = nan(cfg.T,1);
    voltViol = zeros(cfg.T,1);
    slackP = nan(cfg.T,1);
    reverseSubstation = false(cfg.T,1);

    ratings = mpc0.branch(:,RATE_A);
    positiveRated = ratings > 0 & isfinite(ratings);
    hasThermalRatings = any(positiveRated);
    thermalViolCount = 0;
    maxLoadingPct = NaN;

    if hasThermalRatings
        allLoading = [];
    end

    for t=1:cfg.T
        row=sched(t);

        openLines=row.OpenLines;
        if isempty(openLines)
            openLines=determineDefaultOpenLines(mpc0,nOpen);
        end
        capSteps=row.CapSteps;
        if isempty(capSteps)
            capSteps=zeros(1,numel(sys.capBuses));
        end

        radialOK(t)=isRadialAfterOpening(mpc0,openLines);

        m0=prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys);

        % Apply the stored inverter action exactly as in the final replay,
        % including the acceptance safeguard.
        if row.HasDroop
            r0=runpf(m0,mpopt);
            if ~r0.success
                continue;
            end
            V0=abs(r0.bus(:,8)).';
            score0=objectiveVoltage(V0,cfg);

            qmax=computeQmax(Ppv(:,t),mpc0.baseMVA,cfg,sys);
            m1=m0;
            for k=1:numel(sys.pvBuses)
                b=sys.pvBuses(k);
                qpu=pvControl(V0(b),row.Vsp(k),row.Vband(k),qmax(k));
                m1.bus(b,4)=m1.bus(b,4)-qpu*mpc0.baseMVA;
            end

            r1=runpf(m1,mpopt);
            if r1.success
                V1=abs(r1.bus(:,8)).';
                score1=objectiveVoltage(V1,cfg);
                if score1 <= score0
                    r=r1;
                else
                    r=r0;
                end
            else
                r=r0;
            end
        else
            r=runpf(m0,mpopt);
        end

        if ~r.success
            continue;
        end

        pfConverged(t)=true;
        V=abs(r.bus(:,8));
        minV(t)=min(V);
        maxV(t)=max(V);
        voltViol(t)=sum(V<cfg.Vmin | V>cfg.Vmax);

        % Sum active generation at slack/reference buses.
        slackBusIdx=find(r.bus(:,2)==3); % BUS_TYPE=REF
        pSlack=0;
        for sb=slackBusIdx(:).'
            gidx=find(r.gen(:,GEN_BUS)==r.bus(sb,1));
            pSlack=pSlack+sum(r.gen(gidx,PG));
        end
        slackP(t)=pSlack;
        reverseSubstation(t)=pSlack < -1e-9;

        if hasThermalRatings
            Sfrom=sqrt(r.branch(:,PF).^2+r.branch(:,QF).^2);
            Sto=sqrt(r.branch(:,PT).^2+r.branch(:,QT).^2);
            Smax=max(Sfrom,Sto);
            loading=nan(nbr,1);
            loading(positiveRated)=100*Smax(positiveRated)./ratings(positiveRated);
            allLoading=[allLoading; loading(positiveRated)]; %#ok<AGROW>
            thermalViolCount=thermalViolCount+sum(loading(positiveRated)>100+1e-9);
        end
    end

    if hasThermalRatings && ~isempty(allLoading)
        maxLoadingPct=max(allLoading);
    end

    summary=table();
    summary.System=string(sys.name);
    summary.PowerFlowConvergedSteps=sum(pfConverged);
    summary.RadialEndpointSteps=sum(radialOK);
    summary.TotalVoltageViolations=sum(voltViol);
    summary.MinVoltage=min(minV,[],'omitnan');
    summary.MaxVoltage=max(maxV,[],'omitnan');
    summary.PositiveRateABranches=sum(positiveRated);
    summary.ThermalViolations=thermalViolCount;
    summary.MaxRatedBranchLoadingPercent=maxLoadingPct;
    summary.SubstationReverseFlowSteps=sum(reverseSubstation);
    summary.MinSlackActivePowerMW=min(slackP,[],'omitnan');
    summary.MaxSlackActivePowerMW=max(slackP,[],'omitnan');

    R.summary=summary;
    R.pfConverged=pfConverged;
    R.radialOK=radialOK;
    R.voltViol=voltViol;
    R.slackP=slackP;
    R.reverseSubstation=reverseSubstation;
    R.hasThermalRatings=hasThermalRatings;
end


function writeReport(fid,R)
    S=R.summary;
    fprintf(fid,'%s\n',char(S.System));
    fprintf(fid,'%s\n',repmat('-',1,strlength(S.System)));
    fprintf(fid,'AC power-flow convergence: %d/96 steps\n',S.PowerFlowConvergedSteps);
    fprintf(fid,'Connected radial endpoint topology: %d/96 steps\n',S.RadialEndpointSteps);
    fprintf(fid,'Residual voltage-limit bus-time samples: %d\n',S.TotalVoltageViolations);
    fprintf(fid,'Voltage range: %.6f to %.6f p.u.\n',S.MinVoltage,S.MaxVoltage);
    fprintf(fid,'Branches with positive RATE_A data: %d\n',S.PositiveRateABranches);
    if S.PositiveRateABranches>0
        fprintf(fid,'Ex-post rated-branch thermal exceedances: %d\n',S.ThermalViolations);
        fprintf(fid,'Maximum loading among positively rated branches: %.3f %%\n',S.MaxRatedBranchLoadingPercent);
    else
        fprintf(fid,'Thermal diagnostic: NOT AVAILABLE because no positive RATE_A values are present.\n');
    end
    fprintf(fid,'Substation reverse-flow steps (slack P < 0): %d\n',S.SubstationReverseFlowSteps);
    fprintf(fid,'Slack active-power range: %.6f to %.6f MW\n\n', ...
        S.MinSlackActivePowerMW,S.MaxSlackActivePowerMW);
end


function m=prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys)
    m=mpc0;
    m.branch(:,11)=1;
    m.branch(openLines,11)=0;

    P0=m.bus(:,3); Q0=m.bus(:,4);
    m.bus(:,3)=P0.*loadMult(:,t);
    m.bus(:,4)=Q0.*loadMult(:,t);

    for k=1:numel(sys.pvBuses)
        m.bus(sys.pvBuses(k),3)=m.bus(sys.pvBuses(k),3)-Ppv(k,t);
    end
    for k=1:min(numel(sys.capBuses),numel(capSteps))
        m.bus(sys.capBuses(k),4)=m.bus(sys.capBuses(k),4)-capSteps(k)*sys.capStepMVAr;
    end
end


function qmax=computeQmax(PactMW,baseMVA,cfg,sys)
    Smax=sys.pvRatedMW/cfg.PFmin;
    qCircle=sqrt(max(Smax.^2-PactMW.^2,0));
    qPF=PactMW.*tan(acos(cfg.PFmin));
    qmax=min(qCircle,qPF)./baseMVA;
end


function q=pvControl(V,Vsp,Vband,Qmax)
    if Qmax<=0 || Vband<=0
        q=0; return;
    end
    q=-(Qmax/Vband)*(V-Vsp);
    if V-Vsp<=-Vband
        q=Qmax;
    elseif V-Vsp>=Vband
        q=-Qmax;
    end
    q=min(Qmax,max(-Qmax,q));
end


function J=objectiveVoltage(V,cfg)
    under=max(cfg.Vmin-V,0);
    over=max(V-cfg.Vmax,0);
    J=sum((V-1).^2)+1e5*sum(under.^2+over.^2);
end


function ok=isRadialAfterOpening(mpc,openLines)
    nb=size(mpc.bus,1);
    nbr=size(mpc.branch,1);
    status=ones(nbr,1);
    status(openLines)=0;
    closed=find(status==1);

    if numel(closed)~=nb-1
        ok=false; return;
    end

    adj=cell(nb,1);
    for e=closed(:).'
        u=mpc.branch(e,1);
        v=mpc.branch(e,2);
        adj{u}(end+1)=v;
        adj{v}(end+1)=u;
    end

    visited=false(nb,1);
    q=1; visited(1)=true; head=1;
    while head<=numel(q)
        n=q(head); head=head+1;
        for nei=adj{n}
            if ~visited(nei)
                visited(nei)=true;
                q(end+1)=nei;
            end
        end
    end
    ok=all(visited);
end


function openLines=determineDefaultOpenLines(mpc,nOpen)
    existingOpen=find(mpc.branch(:,11)==0).';
    if numel(existingOpen)==nOpen
        openLines=existingOpen;
    else
        openLines=(size(mpc.branch,1)-nOpen+1):size(mpc.branch,1);
    end
end


function rootDir=findProjectRoot(startDir)
    d=startDir;
    while true
        if exist(fullfile(d,'ALEXIS_v4_short_names_bundle'),'dir') || ...
           exist(fullfile(d,'hybrid_dnr_cap_droop_33bus_final_v2_safe.m'),'file')
            rootDir=d; return;
        end
        parent=fileparts(d);
        if strcmp(parent,d)
            rootDir=startDir; return;
        end
        d=parent;
    end
end


function f=findUniqueRecursive(rootDir,fileName,preferredToken)
    D=dir(fullfile(rootDir,'**',fileName));
    if isempty(D)
        error('Could not find %s under %s',fileName,rootDir);
    end
    paths=strings(numel(D),1);
    for k=1:numel(D)
        paths(k)=string(fullfile(D(k).folder,D(k).name));
    end
    if ~isempty(preferredToken)
        hit=contains(paths,preferredToken,'IgnoreCase',true);
        if any(hit), paths=paths(hit); end
    end
    hitV4=contains(paths,'ALEXIS_v4_short_names_bundle','IgnoreCase',true);
    if any(hitV4), paths=paths(hitV4); end
    f=char(paths(1));
end
