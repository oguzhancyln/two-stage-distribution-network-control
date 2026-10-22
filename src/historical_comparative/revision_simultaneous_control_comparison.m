%% revision_simultaneous_control_comparison.m
% Reviewer-revision comparison: simultaneous DNR + capacitor + inverter
% optimization during the PV-active intervals.
%
% This script compares:
%   (i)  the submitted/revised two-stage strategy, which retains the latest
%        topology and capacitor state while adjusting inverter parameters, and
%   (ii) simultaneous control, which is allowed to change topology, capacitor
%        steps, and inverter parameters together at every PV-active time step.
%
% The simultaneous case uses the SAME voltage objective and the SAME dynamic
% inverter capability model. A seed wolf corresponding exactly to the selected
% two-stage state is included at every time step, so the comparison does not
% start from an inferior operating point.
%
% To keep the additional reviewer test computationally practical, the default
% simultaneous-search budget is 30 agents x 30 iterations per PV-active step.
% Change the two settings below if a larger search budget is desired.
%
% Outputs:
%   results_revision_simultaneous/
%       simultaneous_summary_33bus.csv
%       simultaneous_summary_69bus.csv
%       simultaneous_summary_combined.csv
%       simultaneous_report.txt
%       checkpoint_33bus_tXX.mat / checkpoint_69bus_tXX.mat
%
% MATPOWER and case69modified must be on the MATLAB path.

clear; clc;
rng(42);

cfg.T = 96;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambdaUnder = 1e5;
cfg.lambdaOver = 1e5;
cfg.PFmin = 0.90;

cfg.VspLB = 0.95;
cfg.VspUB = 1.05;
cfg.VbandLB = 0.01;
cfg.VbandUB = 0.05;

cfg.nAgents = 30;
cfg.nIter = 30;
cfg.verbose = true;

rootDir = findProjectRoot(pwd);
outDir = fullfile(rootDir,'results_revision_simultaneous');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('Project root: %s\n',rootDir);
fprintf('Simultaneous search budget: %d agents x %d iterations per PV-active step\n\n', ...
    cfg.nAgents,cfg.nIter);

% Final V4 result files provide the two-stage schedules.
f33 = findUniqueRecursive(rootDir,'ablation_results_33bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');
f69 = findUniqueRecursive(rootDir,'ablation_results_69bus.mat', ...
    'results_ablation_33_69_v4_consistent_qmax');

A33 = load(f33);
A69 = load(f69);

sys33 = A33.sys;
sys69 = A69.sys;

% Enforce final paper settings.
sys33.caseName = 'case33bw';
sys33.pvBuses = [21 24 32 17 18];
sys33.capBuses = [18 25 33];
sys33.capStepMVAr = 60/1000;
sys33.capMaxSteps = 5;
sys33.pvRatedMW = 0.20;
sys33.name = '33bus';

sys69.caseName = 'case69modified';
sys69.pvBuses = [15 18 27 35 46 60];
sys69.capBuses = [7 23 30 39 55];
sys69.capStepMVAr = 60/1000;
sys69.capMaxSteps = 5;
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

res33 = runSystem(sys33,A33.result.schedule,load33,pv33,cfg,outDir);
res69 = runSystem(sys69,A69.result.schedule,load69,pv69,cfg,outDir);

writetable(res33.summary,fullfile(outDir,'simultaneous_summary_33bus.csv'));
writetable(res69.summary,fullfile(outDir,'simultaneous_summary_69bus.csv'));
combined = [res33.summary; res69.summary];
writetable(combined,fullfile(outDir,'simultaneous_summary_combined.csv'));

fid = fopen(fullfile(outDir,'simultaneous_report.txt'),'w');
fprintf(fid,'Simultaneous-control comparison\n');
fprintf(fid,'===============================\n\n');
fprintf(fid,'Search budget: %d agents x %d iterations at each PV-active time step.\n', ...
    cfg.nAgents,cfg.nIter);
fprintf(fid,'The first wolf at each step is seeded with the selected two-stage state.\n\n');
writeReport(fid,res33.summary,'33-bus');
writeReport(fid,res69.summary,'69-bus');
fclose(fid);

disp(res33.summary);
disp(res69.summary);

fprintf('\nDone. Please send:\n');
fprintf('  %s\n',fullfile(outDir,'simultaneous_report.txt'));
fprintf('  %s\n',fullfile(outDir,'simultaneous_summary_combined.csv'));


%% ========================================================================
function out = runSystem(sys,sched,loadMult,Ppv,cfg,outDir)

    mpc0 = loadcase(sys.caseName);
    nb = size(mpc0.bus,1);
    nbr = size(mpc0.branch,1);
    nOpen = nbr - (nb-1);
    nCaps = numel(sys.capBuses);
    nPV = numel(sys.pvBuses);

    % First reconstruct the selected two-stage result using the same model.
    Vtwo = nan(cfg.T,nb);
    Vsim = nan(cfg.T,nb);

    selectedOpen = cell(cfg.T,1);
    selectedCaps = cell(cfg.T,1);
    simOpen = cell(cfg.T,1);
    simCaps = cell(cfg.T,1);

    timeRuntime = zeros(cfg.T,1);

    for t = 1:cfg.T
        row = sched(t);

        baseOpen = row.OpenLines;
        if isempty(baseOpen)
            baseOpen = determineDefaultOpenLines(mpc0,nOpen);
        end
        baseCaps = row.CapSteps;
        if isempty(baseCaps)
            baseCaps = zeros(1,nCaps);
        end

        selectedOpen{t} = sort(baseOpen);
        selectedCaps{t} = baseCaps;

        % Reconstruct selected two-stage state.
        if row.HasDroop
            [Vtwo(t,:),~,~,~] = simulateState(mpc0,baseOpen,baseCaps, ...
                true,row.Vsp,row.Vband,t,loadMult,Ppv,sys,cfg,true);
        else
            [Vtwo(t,:),~,~,~] = simulateState(mpc0,baseOpen,baseCaps, ...
                false,[],[],t,loadMult,Ppv,sys,cfg,true);
        end

        if ~row.HasDroop
            % Outside the PV-active intervals, use the same network-control
            % schedule; the comparison concerns simultaneous action only where
            % the two-stage method otherwise adjusts inverter support.
            Vsim(t,:) = Vtwo(t,:);
            simOpen{t} = sort(baseOpen);
            simCaps{t} = baseCaps;
            continue;
        end

        fprintf('\n[%s] simultaneous optimization at t=%d\n',sys.name,t);
        ticStep = tic;

        [bestX,bestScore] = optimizeSimultaneous( ...
            mpc0,baseOpen,baseCaps,row.Vsp,row.Vband, ...
            t,loadMult,Ppv,sys,cfg,nOpen);

        [openBest,capsBest,VspBest,VbandBest,valid] = ...
            decodeCandidate(bestX,mpc0,nOpen,sys,cfg);
        if ~valid
            error('%s t=%d best simultaneous candidate decoded invalid.',sys.name,t);
        end

        [Vbest,~,~,~] = simulateState(mpc0,openBest,capsBest, ...
            true,VspBest,VbandBest,t,loadMult,Ppv,sys,cfg,false);

        Vsim(t,:) = Vbest;
        simOpen{t} = openBest;
        simCaps{t} = capsBest;
        timeRuntime(t) = toc(ticStep);

        save(fullfile(outDir,sprintf('checkpoint_%s_t%02d.mat',sys.name,t)), ...
            't','bestX','bestScore','openBest','capsBest','VspBest','VbandBest','Vbest','timeRuntime');

        fprintf('[%s] t=%d complete | score=%.6g | elapsed=%.1fs | open=%s | caps=%s\n', ...
            sys.name,t,bestScore,timeRuntime(t),mat2str(openBest),mat2str(capsBest));
    end

    metTwo = computeMetrics(Vtwo,cfg);
    metSim = computeMetrics(Vsim,cfg);

    [lineTwo,capTwo] = countControlChanges(selectedOpen,selectedCaps);
    [lineSim,capSim] = countControlChanges(simOpen,simCaps);

    S = table();
    S.System = [string(sys.name); string(sys.name)];
    S.Case = ["Two-stage strategy"; "Simultaneous control"];
    S.MinV = [metTwo.MinV; metSim.MinV];
    S.MaxV = [metTwo.MaxV; metSim.MaxV];
    S.JV = [metTwo.TotalVoltageDeviation; metSim.TotalVoltageDeviation];
    S.UnderVoltageCount = [metTwo.UnderVoltageCount; metSim.UnderVoltageCount];
    S.OverVoltageCount = [metTwo.OverVoltageCount; metSim.OverVoltageCount];
    S.TotalViolationCount = [metTwo.TotalViolationCount; metSim.TotalViolationCount];
    S.LineStateChanges = [lineTwo; lineSim];
    S.CapStepChanges = [capTwo; capSim];
    S.RuntimeSecondsPVActive = [NaN; sum(timeRuntime)];

    out.summary = S;
    out.Vtwo = Vtwo;
    out.Vsim = Vsim;
    out.selectedOpen = selectedOpen;
    out.selectedCaps = selectedCaps;
    out.simOpen = simOpen;
    out.simCaps = simCaps;
    out.runtime = timeRuntime;
end


function [bestX,bestScore] = optimizeSimultaneous( ...
    mpc0,seedOpen,seedCaps,seedVsp,seedVband, ...
    t,loadMult,Ppv,sys,cfg,nOpen)

    nCaps = numel(sys.capBuses);
    nPV = numel(sys.pvBuses);
    nbr = size(mpc0.branch,1);

    dim = nOpen+nCaps+2*nPV;
    lb = [ones(1,nOpen), zeros(1,nCaps), ...
          cfg.VspLB*ones(1,nPV), cfg.VbandLB*ones(1,nPV)];
    ub = [nbr*ones(1,nOpen), sys.capMaxSteps*ones(1,nCaps), ...
          cfg.VspUB*ones(1,nPV), cfg.VbandUB*ones(1,nPV)];

    % Some saved sys structs may not retain capMaxSteps.
    if ~isfield(sys,'capMaxSteps') || isempty(sys.capMaxSteps)
        sys.capMaxSteps = 5;
        ub = [nbr*ones(1,nOpen), sys.capMaxSteps*ones(1,nCaps), ...
              cfg.VspUB*ones(1,nPV), cfg.VbandUB*ones(1,nPV)];
    end

    X = zeros(cfg.nAgents,dim);

    % Wolf 1 = exact selected two-stage state.
    X(1,:) = [seedOpen(:).',seedCaps(:).',seedVsp(:).',seedVband(:).'];

    % Other wolves start from feasible radial topologies near the selected one.
    for i = 2:cfg.nAgents
        o = seedOpen;
        nExchange = randi([1,3]);
        for q = 1:nExchange
            o = randomBranchExchange(mpc0,o);
        end
        c = randi([0,sys.capMaxSteps],1,nCaps);
        vsp = cfg.VspLB + rand(1,nPV)*(cfg.VspUB-cfg.VspLB);
        vb  = cfg.VbandLB + rand(1,nPV)*(cfg.VbandUB-cfg.VbandLB);
        X(i,:) = [o,c,vsp,vb];
    end

    alphaScore=inf; betaScore=inf; deltaScore=inf;
    alphaX=zeros(1,dim); betaX=alphaX; deltaX=alphaX;

    for it = 1:cfg.nIter

        nFeasible = 0;
        for i = 1:cfg.nAgents
            X(i,:) = reflectBounds(X(i,:),lb,ub);

            [openLines,capSteps,Vsp,Vband,valid] = ...
                decodeCandidate(X(i,:),mpc0,nOpen,sys,cfg);

            if ~valid
                score = 1e12;
            else
                nFeasible = nFeasible+1;
                [~,score,~,~] = simulateState(mpc0,openLines,capSteps, ...
                    true,Vsp,Vband,t,loadMult,Ppv,sys,cfg,false);
                if ~isfinite(score), score=1e12; end
            end

            if score < alphaScore
                deltaScore=betaScore; deltaX=betaX;
                betaScore=alphaScore; betaX=alphaX;
                alphaScore=score; alphaX=X(i,:);
            elseif score < betaScore
                deltaScore=betaScore; deltaX=betaX;
                betaScore=score; betaX=X(i,:);
            elseif score < deltaScore
                deltaScore=score; deltaX=X(i,:);
            end
        end

        a = 2-it*(2/cfg.nIter);
        for i = 1:cfg.nAgents
            for j = 1:dim
                r1=rand; r2=rand; A1=2*a*r1-a; C1=2*r2;
                X1=alphaX(j)-A1*abs(C1*alphaX(j)-X(i,j));

                r1=rand; r2=rand; A2=2*a*r1-a; C2=2*r2;
                X2=betaX(j)-A2*abs(C2*betaX(j)-X(i,j));

                r1=rand; r2=rand; A3=2*a*r1-a; C3=2*r2;
                X3=deltaX(j)-A3*abs(C3*deltaX(j)-X(i,j));

                X(i,j)=(X1+X2+X3)/3;
            end
            X(i,:)=reflectBounds(X(i,:),lb,ub);
        end

        if cfg.verbose
            fprintf('  iter %2d/%2d | best=%12.6g | feasible=%2d/%2d\n', ...
                it,cfg.nIter,alphaScore,nFeasible,cfg.nAgents);
        end
    end

    bestX=alphaX;
    bestScore=alphaScore;
end


function [openLines,capSteps,Vsp,Vband,valid] = decodeCandidate(x,mpc,nOpen,sys,cfg)

    nbr = size(mpc.branch,1);
    nCaps = numel(sys.capBuses);
    nPV = numel(sys.pvBuses);

    rawOpen = round(x(1:nOpen));
    rawOpen = min(nbr,max(1,rawOpen));
    openLines = sort(unique(rawOpen,'stable'));

    capStart = nOpen+1;
    capEnd = nOpen+nCaps;
    capSteps = round(x(capStart:capEnd));
    capSteps = min(sys.capMaxSteps,max(0,capSteps));

    vspStart = capEnd+1;
    vspEnd = capEnd+nPV;
    Vsp = min(cfg.VspUB,max(cfg.VspLB,x(vspStart:vspEnd)));

    vbStart = vspEnd+1;
    vbEnd = vspEnd+nPV;
    Vband = min(cfg.VbandUB,max(cfg.VbandLB,x(vbStart:vbEnd)));

    valid = numel(openLines)==nOpen && isRadialAfterOpening(mpc,openLines);
    if ~valid
        openLines=[];
    end
end


function [Vfinal,score,accepted,rejected] = simulateState( ...
    mpc0,openLines,capSteps,useDroop,Vsp,Vband,t,loadMult,Ppv,sys,cfg,gate)

    accepted=0; rejected=0;
    m0=prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys);
    mpopt=mpoption('verbose',0,'out.all',0);

    if ~useDroop
        r=runpf(m0,mpopt);
        if ~r.success
            Vfinal=nan(1,size(mpc0.bus,1)); score=inf; return;
        end
        Vfinal=abs(r.bus(:,8)).';
        score=objectiveVoltage(Vfinal,cfg);
        return;
    end

    r0=runpf(m0,mpopt);
    if ~r0.success
        Vfinal=nan(1,size(mpc0.bus,1)); score=inf; return;
    end
    V0=abs(r0.bus(:,8)).';
    score0=objectiveVoltage(V0,cfg);

    QmaxStep=computeQmax(Ppv(:,t),mpc0.baseMVA,cfg,sys);
    m1=m0;
    for k=1:numel(sys.pvBuses)
        b=sys.pvBuses(k);
        qpu=pvControl(V0(b),Vsp(k),Vband(k),QmaxStep(k));
        m1.bus(b,4)=m1.bus(b,4)-qpu*mpc0.baseMVA;
    end

    r1=runpf(m1,mpopt);
    if ~r1.success
        Vfinal=V0; score=score0; rejected=1; return;
    end

    V1=abs(r1.bus(:,8)).';
    score1=objectiveVoltage(V1,cfg);

    if gate && score1>score0
        Vfinal=V0; score=score0; rejected=1;
    else
        Vfinal=V1; score=score1; accepted=1;
    end
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


function QmaxStep=computeQmax(PactMW,baseMVA,cfg,sys)
    SmaxMW=sys.pvRatedMW/cfg.PFmin;
    qCircle=sqrt(max(SmaxMW.^2-PactMW.^2,0));
    qPF=PactMW.*tan(acos(cfg.PFmin));
    qmax=min(qCircle,qPF);
    QmaxStep=qmax/baseMVA;
end


function q=pvControl(V,Vsp,Vband,Qmax)
    if Qmax<=0 || Vband<=0
        q=0; return;
    end
    K=Qmax/Vband;
    dV=V-Vsp;
    q=-K*dV;
    if dV<=-Vband
        q=Qmax;
    elseif dV>=Vband
        q=-Qmax;
    end
    q=min(Qmax,max(-Qmax,q));
end


function J=objectiveVoltage(V,cfg)
    if any(isnan(V)) || any(isinf(V))
        J=inf; return;
    end
    under=max(cfg.Vmin-V,0);
    over=max(V-cfg.Vmax,0);
    J=sum((V-1).^2)+cfg.lambdaUnder*sum(under.^2)+cfg.lambdaOver*sum(over.^2);
end


function y=reflectBounds(x,lb,ub)
    y=x;
    for j=1:numel(x)
        lo=lb(j); hi=ub(j);
        if hi<=lo
            y(j)=lo; continue;
        end
        d=hi-lo;
        y(j)=lo+abs(mod(y(j)-hi,2*d)-d);
        y(j)=min(hi,max(lo,y(j)));
    end
end


function newOpen=randomBranchExchange(mpc,openLines)
    % Starting from a radial topology, close one open branch to create a cycle,
    % then open one branch on the pre-existing path between its endpoints.
    openLines=sort(openLines);
    nbr=size(mpc.branch,1);

    tie=openLines(randi(numel(openLines)));
    u=mpc.branch(tie,1);
    v=mpc.branch(tie,2);

    closed=setdiff(1:nbr,openLines);
    [pathBranches,ok]=findPathBranches(mpc,closed,u,v);
    if ~ok || isempty(pathBranches)
        newOpen=openLines;
        return;
    end

    branchToOpen=pathBranches(randi(numel(pathBranches)));
    newOpen=sort([setdiff(openLines,tie),branchToOpen]);

    if ~isRadialAfterOpening(mpc,newOpen)
        newOpen=openLines;
    end
end


function [pathBranches,ok]=findPathBranches(mpc,closed,u,v)
    nb=size(mpc.bus,1);
    adj=cell(nb,1);
    edgeIdx=cell(nb,1);

    for e=closed
        a=mpc.branch(e,1);
        b=mpc.branch(e,2);
        adj{a}(end+1)=b;
        edgeIdx{a}(end+1)=e;
        adj{b}(end+1)=a;
        edgeIdx{b}(end+1)=e;
    end

    parent=zeros(nb,1);
    parentEdge=zeros(nb,1);
    visited=false(nb,1);
    q=u; visited(u)=true; head=1;

    while head<=numel(q)
        n=q(head); head=head+1;
        if n==v, break; end
        for kk=1:numel(adj{n})
            nei=adj{n}(kk);
            if ~visited(nei)
                visited(nei)=true;
                parent(nei)=n;
                parentEdge(nei)=edgeIdx{n}(kk);
                q(end+1)=nei;
            end
        end
    end

    if ~visited(v)
        pathBranches=[]; ok=false; return;
    end

    pathBranches=[];
    cur=v;
    while cur~=u
        pathBranches(end+1)=parentEdge(cur);
        cur=parent(cur);
    end
    ok=true;
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
        u=mpc.branch(e,1); v=mpc.branch(e,2);
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


function [lineChanges,capChanges]=countControlChanges(openCell,capCell)
    lineChanges=0;
    capChanges=0;
    for t=2:numel(openCell)
        if isempty(openCell{t-1}) || isempty(openCell{t}), continue; end
        lineChanges=lineChanges+numel(setxor(openCell{t-1},openCell{t}));
        if ~isempty(capCell{t-1}) && ~isempty(capCell{t})
            capChanges=capChanges+sum(abs(capCell{t}-capCell{t-1}));
        end
    end
end


function met=computeMetrics(Vmat,cfg)
    V=Vmat(:);
    V=V(~isnan(V));
    met.MinV=min(V);
    met.MaxV=max(V);
    met.TotalVoltageDeviation=sum((V-1).^2);
    met.UnderVoltageCount=sum(V<cfg.Vmin);
    met.OverVoltageCount=sum(V>cfg.Vmax);
    met.TotalViolationCount=met.UnderVoltageCount+met.OverVoltageCount;
end


function openLines=determineDefaultOpenLines(mpc,nOpen)
    existingOpen=find(mpc.branch(:,11)==0).';
    if numel(existingOpen)==nOpen
        openLines=existingOpen;
    else
        openLines=(size(mpc.branch,1)-nOpen+1):size(mpc.branch,1);
    end
end


function writeReport(fid,S,name)
    fprintf(fid,'%s\n',name);
    fprintf(fid,'%s\n',repmat('-',1,length(name)));
    for i=1:height(S)
        fprintf(fid,'%-22s MinV=%8.5f MaxV=%8.5f JV=%9.4f Viol=%4d LineChanges=%4d CapChanges=%4d Runtime=%9.1f\n', ...
            char(S.Case(i)),S.MinV(i),S.MaxV(i),S.JV(i), ...
            S.TotalViolationCount(i),S.LineStateChanges(i),S.CapStepChanges(i), ...
            S.RuntimeSecondsPVActive(i));
    end
    fprintf(fid,'\n');
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
