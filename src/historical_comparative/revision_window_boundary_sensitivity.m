%% revision_window_boundary_sensitivity.m
% Reviewer-revision experiment: control-window boundary sensitivity.
%
% Purpose
% -------
% Re-optimize the complete two-stage workflow for three schedules whose
% internal boundaries are shifted by +/-2 time steps (30 min):
%
%   Sminus2 : network [1-22], [38-58], [71-96]
%             inverter [23-37], [59-70]
%   S0      : network [1-24], [40-60], [73-96]
%             inverter [25-39], [61-72]
%   Splus2  : network [1-26], [42-62], [75-96]
%             inverter [27-41], [63-74]
%
% IMPORTANT INTERPRETATION
% ------------------------
% This is an equal-budget boundary-sensitivity experiment, not a rerun of the
% headline 50x100 / 40x50 solution. To keep the reviewer test computationally
% manageable, every schedule uses exactly the same reduced search budget:
%   network control : 30 wolves x 60 iterations
%   inverter tuning : 30 wolves x 30 iterations
%
% The same random seed is reset at the start of each schedule within a feeder
% so that the comparison isolates schedule-boundary effects as far as possible.
% Do NOT compare the absolute S0 values from this experiment directly with the
% headline selected-schedule values; compare Sminus2, S0, and Splus2 with each
% other under this common reduced budget.
%
% Final paper-consistent modeling choices
% ---------------------------------------
% - MATPOWER AC power flow
% - connected radial endpoint topology
% - continuous GWO positions decoded to mixed variables
% - invalid decoded topology -> large penalty; no repair after GWO update
% - capacitor steps: integer 0..5, 60 kVAr/step
% - PFmin = 0.90
% - Qmax = min(sqrt(S^2-P^2), P*tan(acos(PFmin)))
% - one-shot inverter Q action followed by a second AC power flow
% - inverter action accepted only if Jpost <= Jpre
% - voltage limits 0.95--1.05 p.u. are soft objective penalties
%
% Outputs
% -------
% results_revision_window_sensitivity/
%   window_sensitivity_summary.csv
%   window_sensitivity_report.txt
%   sensitivity_<system>_<schedule>.mat
%   sensitivity_schedule_<system>_<schedule>.csv
%
% Recommended run
% ---------------
% Put this file in the ALEXIS_PAPER project root and run:
%   run('revision_window_boundary_sensitivity.m')
%
% Expected wall-clock time on the previously benchmarked laptop is several
% hours. Use an uninterrupted session if possible.

clear; clc; close all;
warning('off','all');

%% -------------------------- SETTINGS -----------------------------------
cfg.T = 96;
cfg.dtHours = 0.25;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambda = 1e5;
cfg.PFmin = 0.90;
cfg.maxCapStep = 5;
cfg.capStepMVAr = 60/1000;

% Equal reduced budget for all three schedules.
cfg.netAgents = 30;
cfg.netIter = 60;
cfg.invAgents = 30;
cfg.invIter = 30;
cfg.verboseNet = true;
cfg.verboseInv = false;

% Same seed across Sminus2/S0/Splus2 within each feeder.
cfg.seed33 = 77001;
cfg.seed69 = 78001;

cfg.schedules(1).name = 'Sminus2';
cfg.schedules(1).network = [1 22; 38 58; 71 96];
cfg.schedules(1).inverter = [23 37; 59 70];

cfg.schedules(2).name = 'S0';
cfg.schedules(2).network = [1 24; 40 60; 73 96];
cfg.schedules(2).inverter = [25 39; 61 72];

cfg.schedules(3).name = 'Splus2';
cfg.schedules(3).network = [1 26; 42 62; 75 96];
cfg.schedules(3).inverter = [27 41; 63 74];

rootDir = findProjectRoot(pwd);
outDir = fullfile(rootDir,'results_revision_window_sensitivity');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\n============================================================\n');
fprintf(' WINDOW-BOUNDARY SENSITIVITY | EQUAL REDUCED BUDGET\n');
fprintf('============================================================\n');
fprintf('Network GWO : %d wolves x %d iterations\n',cfg.netAgents,cfg.netIter);
fprintf('Inverter GWO: %d wolves x %d iterations\n',cfg.invAgents,cfg.invIter);
fprintf('PFmin       : %.2f\n',cfg.PFmin);
fprintf('Project root: %s\n',rootDir);
fprintf(['Interpretation: compare Sminus2/S0/Splus2 only with each other; ' ...
         'do not compare their absolute values directly with headline results.\n\n']);

systems = makeSystemConfigs(rootDir,cfg);
summaryRows = {};

for sidx = 1:numel(systems)
    sys = systems(sidx);
    fprintf('\n================ %s sensitivity ================\n',sys.name);

    [mpc0,Lprof,PVprof,defaultOpen] = loadSystemAndProfiles(sys,cfg);
    fprintf('Loaded %s: %d buses, %d branches, default open=%s\n', ...
        sys.caseName,size(mpc0.bus,1),size(mpc0.branch,1),mat2str(defaultOpen));
    fprintf('Profiles: load range=[%.3f, %.3f], PV max=%.4f MW/site\n', ...
        min(Lprof(:)),max(Lprof(:)),max(PVprof(:)));

    for q = 1:numel(cfg.schedules)
        schedDef = cfg.schedules(q);
        if strcmpi(sys.name,'33bus'), seed=cfg.seed33; else, seed=cfg.seed69; end
        rng(seed,'twister');

        fprintf('\n------------------------------------------------------------\n');
        fprintf('%s | %s | seed=%d\n',sys.name,schedDef.name,seed);
        fprintf('Network : [%d-%d], [%d-%d], [%d-%d]\n', ...
            schedDef.network(1,1),schedDef.network(1,2), ...
            schedDef.network(2,1),schedDef.network(2,2), ...
            schedDef.network(3,1),schedDef.network(3,2));
        fprintf('Inverter: [%d-%d], [%d-%d]\n', ...
            schedDef.inverter(1,1),schedDef.inverter(1,2), ...
            schedDef.inverter(2,1),schedDef.inverter(2,2));

        runTimer = tic;
        result = optimizeSchedule(mpc0,Lprof,PVprof,defaultOpen,sys,cfg,schedDef,outDir);
        runtime = toc(runTimer);

        M = computeMetrics(result.V,cfg);
        [lineChanges,capChanges] = countDiscreteChanges(defaultOpen, ...
            zeros(1,numel(sys.capBuses)),result.openByT,result.capsByT);
        activeLoss_kWh = sum(result.PlossMW)*cfg.dtHours*1000;

        summaryRows(end+1,:) = {sys.name,schedDef.name,seed, ...
            windowsToString(schedDef.network),windowsToString(schedDef.inverter), ...
            M.TotalVoltageDeviation,M.TotalViolationCount,M.UnderVoltageCount, ...
            M.OverVoltageCount,M.MinV,M.MaxV,activeLoss_kWh, ...
            lineChanges,capChanges,result.acceptedCount,result.rejectedCount,runtime}; %#ok<SAGROW>

        save(fullfile(outDir,sprintf('sensitivity_%s_%s.mat',sys.name,schedDef.name)), ...
            'sys','cfg','schedDef','seed','result','M','lineChanges','capChanges', ...
            'activeLoss_kWh','runtime');
        writeScheduleCSV(result,sys,schedDef,outDir);

        fprintf('%s %s DONE | JV=%.6f | viol=%d | min=%.5f | max=%.5f | loss=%.1f kWh\n', ...
            sys.name,schedDef.name,M.TotalVoltageDeviation,M.TotalViolationCount, ...
            M.MinV,M.MaxV,activeLoss_kWh);
        fprintf('Discrete changes: line=%d | cap steps=%d | inverter accepted/rejected=%d/%d\n', ...
            lineChanges,capChanges,result.acceptedCount,result.rejectedCount);
        fprintf('Runtime: %.1f s (%.1f min)\n',runtime,runtime/60);

        % Write cumulative summary after every completed schedule.
        Tpartial = rowsToTable(summaryRows);
        writetable(Tpartial,fullfile(outDir,'window_sensitivity_summary_partial.csv'));
    end
end

T = rowsToTable(summaryRows);
writetable(T,fullfile(outDir,'window_sensitivity_summary.csv'));
writeReport(T,cfg,fullfile(outDir,'window_sensitivity_report.txt'));

fprintf('\n============================================================\n');
fprintf('DONE. Send these two files for manuscript integration:\n');
fprintf('  %s\n',fullfile(outDir,'window_sensitivity_summary.csv'));
fprintf('  %s\n',fullfile(outDir,'window_sensitivity_report.txt'));
fprintf('============================================================\n');


%% ========================== MAIN WORKFLOW ===============================
function result = optimizeSchedule(mpc0,Lprof,PVprof,defaultOpen,sys,cfg,schedDef,outDir)
    T = cfg.T;
    nb = size(mpc0.bus,1);
    nCaps = numel(sys.capBuses);
    nPV = numel(sys.pvBuses);

    openByT = repmat({sort(defaultOpen)},T,1);
    capsByT = repmat({zeros(1,nCaps)},T,1);
    VspByT = repmat({nan(1,nPV)},T,1);
    VbandByT = repmat({nan(1,nPV)},T,1);
    accepted = false(T,1);
    evaluated = false(T,1);

    % ---------------- Stage 1: all three network-control windows ----------
    for w = 1:3
        tr = schedDef.network(w,1):schedDef.network(w,2);
        label = sprintf('%s %s NET W%d t=%d..%d',sys.name,schedDef.name,w,tr(1),tr(end));
        [bestOpen,bestCaps,bestObj,conv] = gwoNetwork( ...
            mpc0,Lprof,PVprof,sys,cfg,tr,label);

        for t = tr
            openByT{t} = sort(bestOpen);
            capsByT{t} = bestCaps;
        end

        % Carry W1 and W2 states into their following inverter intervals.
        if w <= 2
            dr = schedDef.inverter(w,1):schedDef.inverter(w,2);
            for t = dr
                openByT{t} = sort(bestOpen);
                capsByT{t} = bestCaps;
            end
        end

        checkpoint = struct('window',w,'tRange',tr,'bestOpen',bestOpen, ...
            'bestCaps',bestCaps,'bestObj',bestObj,'convergence',conv);
        save(fullfile(outDir,sprintf('checkpoint_%s_%s_netW%d.mat', ...
            sys.name,schedDef.name,w)),'checkpoint');

        fprintf('  W%d selected | obj=%.6g | open=%s | caps=%s\n', ...
            w,bestObj,mat2str(bestOpen),mat2str(bestCaps));
    end

    % -------- Evaluate retained network-control states over 96 steps ------
    Vnetwork = nan(T,nb);
    for t=1:T
        [Vnetwork(t,:),~,ok] = simulateNoInverter( ...
            mpc0,Lprof,PVprof,sys,cfg,t,openByT{t},capsByT{t});
        if ~ok
            error('%s %s: AC power flow failed at network-only t=%d.', ...
                sys.name,schedDef.name,t);
        end
    end

    % ---------------- Stage 2: inverter tuning ----------------------------
    Vfinal = Vnetwork;
    PlossMW = nan(T,1);
    acceptedCount = 0;
    rejectedCount = 0;

    for t=1:T
        if isInWindows(t,schedDef.inverter)
            evaluated(t)=true;
            label = sprintf('%s %s INV t=%d',sys.name,schedDef.name,t);
            [bestEta,~,~,~] = gwoInverterOneStep( ...
                mpc0,Lprof,PVprof,sys,cfg,t,openByT{t},capsByT{t},label);

            nPV = numel(sys.pvBuses);
            VspByT{t}=bestEta(1:nPV);
            VbandByT{t}=bestEta(nPV+1:end);

            [Vpre,Vpost,Jpre,Jpost,lossPre,lossPost,okPost] = ...
                simulateWithInverter(mpc0,Lprof,PVprof,sys,cfg,t, ...
                openByT{t},capsByT{t},bestEta);

            if okPost && Jpost <= Jpre
                Vfinal(t,:)=Vpost;
                PlossMW(t)=lossPost;
                accepted(t)=true;
                acceptedCount=acceptedCount+1;
            else
                Vfinal(t,:)=Vpre;
                PlossMW(t)=lossPre;
                rejectedCount=rejectedCount+1;
            end

            fprintf('  inverter t=%2d | %s | Jpre=%.6g Jpost=%.6g\n', ...
                t,ternary(accepted(t),'ACCEPT','REJECT'),Jpre,Jpost);
        else
            [V,loss,ok] = simulateNoInverter( ...
                mpc0,Lprof,PVprof,sys,cfg,t,openByT{t},capsByT{t});
            if ~ok, error('Final PF failed at t=%d.',t); end
            Vfinal(t,:)=V;
            PlossMW(t)=loss;
        end
    end

    result = struct();
    result.V = Vfinal;
    result.Vnetwork = Vnetwork;
    result.PlossMW = PlossMW;
    result.openByT = openByT;
    result.capsByT = capsByT;
    result.VspByT = VspByT;
    result.VbandByT = VbandByT;
    result.accepted = accepted;
    result.evaluated = evaluated;
    result.acceptedCount = acceptedCount;
    result.rejectedCount = rejectedCount;
end


%% ============================ NETWORK GWO ===============================
function [bestOpen,bestCaps,bestObj,conv] = gwoNetwork(mpc0,Lprof,PVprof,sys,cfg,tRange,label)
    nb = size(mpc0.bus,1);
    nBr = size(mpc0.branch,1);
    nOpen = nBr-(nb-1);
    nCaps = numel(sys.capBuses);
    dim = nOpen+nCaps;

    lb = [ones(1,nOpen),zeros(1,nCaps)];
    ub = [nBr*ones(1,nOpen),cfg.maxCapStep*ones(1,nCaps)];

    X=zeros(cfg.netAgents,dim);
    for a=1:cfg.netAgents
        X(a,1:nOpen)=randomRadialOpen(mpc0,nOpen);
        X(a,nOpen+1:end)=randi([0 cfg.maxCapStep],1,nCaps);
    end

    alpha=zeros(1,dim); beta=alpha; delta=alpha;
    alphaScore=inf; betaScore=inf; deltaScore=inf;
    conv=nan(cfg.netIter,1);

    for it=1:cfg.netIter
        feasible=0;
        for a=1:cfg.netAgents
            X(a,:)=reflectBounds(X(a,:),lb,ub);
            [openLines,caps,valid]=decodeNetworkCandidate( ...
                X(a,:),mpc0,nOpen,nCaps,cfg);
            if ~valid
                score=1e12;
            else
                feasible=feasible+1;
                score=networkFitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps);
            end

            if score<alphaScore
                deltaScore=betaScore; delta=beta;
                betaScore=alphaScore; beta=alpha;
                alphaScore=score; alpha=X(a,:);
            elseif score<betaScore
                deltaScore=betaScore; delta=beta;
                betaScore=score; beta=X(a,:);
            elseif score<deltaScore
                deltaScore=score; delta=X(a,:);
            end
        end

        aGWO=2-it*(2/cfg.netIter);
        X=updateGWO(X,alpha,beta,delta,aGWO,lb,ub);
        conv(it)=alphaScore;

        if cfg.verboseNet && (it==1 || mod(it,10)==0 || it==cfg.netIter)
            fprintf('[%s] iter %2d/%2d | best=%12.6g | feasible=%d/%d\n', ...
                label,it,cfg.netIter,alphaScore,feasible,cfg.netAgents);
        end
    end

    [bestOpen,bestCaps,valid]=decodeNetworkCandidate(alpha,mpc0,nOpen,nCaps,cfg);
    if ~valid
        error('%s: best network candidate is invalid.',label);
    end
    bestObj=alphaScore;
end

function score=networkFitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps)
    score=0;
    for t=tRange
        [V,~,ok]=simulateNoInverter(mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps);
        if ~ok
            score=score+1e10;
        else
            score=score+voltageObjective(V,cfg);
        end
    end
end

function [openLines,caps,valid]=decodeNetworkCandidate(x,mpc0,nOpen,nCaps,cfg)
    nBr=size(mpc0.branch,1);
    candidate=round(x(1:nOpen));
    candidate=min(nBr,max(1,candidate));
    openLines=unique(candidate,'stable');
    openLines=sort(openLines);

    capRaw=round(x(nOpen+1:nOpen+nCaps));
    caps=min(cfg.maxCapStep,max(0,capRaw));

    valid=(numel(openLines)==nOpen) && isRadialAfterOpening(mpc0,openLines);
    if ~valid, openLines=[]; end
end

function openLines=randomRadialOpen(mpc0,nOpen)
    nBr=size(mpc0.branch,1);
    for tries=1:20000
        cand=sort(randperm(nBr,nOpen));
        if isRadialAfterOpening(mpc0,cand)
            openLines=cand;
            return;
        end
    end
    error('Could not generate a connected radial opening set.');
end

function ok=isRadialAfterOpening(mpc0,openLines)
    nb=size(mpc0.bus,1);
    nBr=size(mpc0.branch,1);
    status=ones(nBr,1);
    status(openLines)=0;
    closed=find(status==1);
    if numel(closed)~=nb-1
        ok=false; return;
    end

    adj=cell(nb,1);
    for e=closed(:).'
        u=mpc0.branch(e,1); v=mpc0.branch(e,2);
        adj{u}(end+1)=v;
        adj{v}(end+1)=u;
    end
    visited=false(nb,1);
    q=1; visited(1)=true; head=1;
    while head<=numel(q)
        u=q(head); head=head+1;
        for v=adj{u}
            if ~visited(v)
                visited(v)=true;
                q(end+1)=v; %#ok<AGROW>
            end
        end
    end
    ok=all(visited);
end


%% ========================== INVERTER GWO ================================
function [bestEta,bestObj,bestV,conv]=gwoInverterOneStep( ...
    mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps,label)

    nPV=numel(sys.pvBuses);
    dim=2*nPV;
    lb=[0.95*ones(1,nPV),0.01*ones(1,nPV)];
    ub=[1.05*ones(1,nPV),0.05*ones(1,nPV)];

    X=rand(cfg.invAgents,dim).*(ub-lb)+lb;
    alpha=zeros(1,dim); beta=alpha; delta=alpha;
    alphaScore=inf; betaScore=inf; deltaScore=inf;
    bestV=nan(1,size(mpc0.bus,1));
    conv=nan(cfg.invIter,1);

    for it=1:cfg.invIter
        for a=1:cfg.invAgents
            X(a,:)=max(min(X(a,:),ub),lb);
            [~,Vpost,~,Jpost,~,~,okPost]=simulateWithInverter( ...
                mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps,X(a,:));
            if ~okPost, score=1e10; else, score=Jpost; end

            if score<alphaScore
                deltaScore=betaScore; delta=beta;
                betaScore=alphaScore; beta=alpha;
                alphaScore=score; alpha=X(a,:); bestV=Vpost;
            elseif score<betaScore
                deltaScore=betaScore; delta=beta;
                betaScore=score; beta=X(a,:);
            elseif score<deltaScore
                deltaScore=score; delta=X(a,:);
            end
        end

        aGWO=2-it*(2/cfg.invIter);
        X=updateGWO(X,alpha,beta,delta,aGWO,lb,ub);
        conv(it)=alphaScore;

        if cfg.verboseInv && (it==1 || mod(it,10)==0 || it==cfg.invIter)
            fprintf('[%s] iter %2d/%2d | best=%12.6g\n', ...
                label,it,cfg.invIter,alphaScore);
        end
    end
    bestEta=alpha;
    bestObj=alphaScore;
end

function Xnew=updateGWO(X,alpha,beta,delta,aGWO,lb,ub)
    Xnew=X;
    for i=1:size(X,1)
        for j=1:size(X,2)
            r1=rand; r2=rand; A1=2*aGWO*r1-aGWO; C1=2*r2;
            X1=alpha(j)-A1*abs(C1*alpha(j)-X(i,j));
            r1=rand; r2=rand; A2=2*aGWO*r1-aGWO; C2=2*r2;
            X2=beta(j)-A2*abs(C2*beta(j)-X(i,j));
            r1=rand; r2=rand; A3=2*aGWO*r1-aGWO; C3=2*r2;
            X3=delta(j)-A3*abs(C3*delta(j)-X(i,j));
            Xnew(i,j)=(X1+X2+X3)/3;
        end
    end
    Xnew=reflectBounds(Xnew,lb,ub);
end

function y=reflectBounds(x,lb,ub)
    y=x;
    if isvector(y), y=reshape(y,1,[]); end
    for i=1:size(y,1)
        for j=1:size(y,2)
            lo=lb(j); hi=ub(j); d=hi-lo;
            if d<=0
                y(i,j)=lo;
            else
                y(i,j)=lo+abs(mod(y(i,j)-hi,2*d)-d);
                y(i,j)=min(hi,max(lo,y(i,j)));
            end
        end
    end
    if size(x,1)==1 && isvector(x), y=y(1,:); end
end


%% ========================== POWER FLOW MODEL ============================
function [V,PlossMW,ok]=simulateNoInverter( ...
    mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps)

    mpc=buildCase(mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps);
    mpopt=mpoption('verbose',0,'out.all',0);
    r=runpf(mpc,mpopt);
    ok=r.success;
    if ~ok
        V=0.5*ones(1,size(mpc0.bus,1));
        PlossMW=NaN;
        return;
    end
    V=abs(r.bus(:,8)).';
    PlossMW=sum(r.branch(:,14)+r.branch(:,16));
end

function [Vpre,Vpost,Jpre,Jpost,lossPre,lossPost,okPost]=simulateWithInverter( ...
    mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps,eta)

    mpc=buildCase(mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps);
    mpopt=mpoption('verbose',0,'out.all',0);
    r0=runpf(mpc,mpopt);
    if ~r0.success
        Vpre=0.5*ones(1,size(mpc0.bus,1));
        Vpost=Vpre; Jpre=1e10; Jpost=1e10;
        lossPre=NaN; lossPost=NaN; okPost=false; return;
    end

    Vpre=abs(r0.bus(:,8)).';
    Jpre=voltageObjective(Vpre,cfg);
    lossPre=sum(r0.branch(:,14)+r0.branch(:,16));

    nPV=numel(sys.pvBuses);
    Vsp=eta(1:nPV);
    Vband=eta(nPV+1:end);
    QmaxStep=computeQmax(PVprof(:,t),mpc0.baseMVA,cfg,sys);

    m1=mpc;
    for k=1:nPV
        b=sys.pvBuses(k);
        qpu=pvControl(Vpre(b),Vsp(k),Vband(k),QmaxStep(k));
        m1.bus(b,4)=m1.bus(b,4)-qpu*mpc0.baseMVA;
    end

    r1=runpf(m1,mpopt);
    okPost=r1.success;
    if ~okPost
        Vpost=Vpre; Jpost=1e10; lossPost=NaN; return;
    end
    Vpost=abs(r1.bus(:,8)).';
    Jpost=voltageObjective(Vpost,cfg);
    lossPost=sum(r1.branch(:,14)+r1.branch(:,16));
end

function mpc=buildCase(mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps)
    mpc=mpc0;
    mpc.branch(:,11)=1;
    mpc.branch(openLines,11)=0;

    P0=mpc.bus(:,3); Q0=mpc.bus(:,4);
    mpc.bus(:,3)=P0.*Lprof(:,t);
    mpc.bus(:,4)=Q0.*Lprof(:,t);

    for k=1:numel(sys.pvBuses)
        mpc.bus(sys.pvBuses(k),3)=mpc.bus(sys.pvBuses(k),3)-PVprof(k,t);
    end

    for c=1:numel(sys.capBuses)
        mpc.bus(sys.capBuses(c),4)=mpc.bus(sys.capBuses(c),4)-caps(c)*cfg.capStepMVAr;
    end
end

function QmaxStep=computeQmax(PactMW,baseMVA,cfg,sys)
    SmaxMW=sys.pvRatedMW./cfg.PFmin;
    qCircle=sqrt(max(SmaxMW.^2-PactMW.^2,0));
    qPF=PactMW.*tan(acos(cfg.PFmin));
    QmaxStep=min(qCircle,qPF)./baseMVA;
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

function J=voltageObjective(V,cfg)
    under=max(cfg.Vmin-V,0);
    over=max(V-cfg.Vmax,0);
    J=sum((V-1).^2)+cfg.lambda*sum(under.^2+over.^2);
end


%% ========================== METRICS / OUTPUT ============================
function M=computeMetrics(V,cfg)
    under=V<cfg.Vmin;
    over=V>cfg.Vmax;
    M.MinV=min(V(:));
    M.MaxV=max(V(:));
    M.MeanV=mean(V(:));
    M.TotalVoltageDeviation=sum((V(:)-1).^2);
    M.UnderVoltageCount=nnz(under);
    M.OverVoltageCount=nnz(over);
    M.TotalViolationCount=nnz(under|over);
end

function [lineChanges,capChanges]=countDiscreteChanges(defaultOpen,defaultCaps,openByT,capsByT)
    lineChanges=numel(setxor(sort(defaultOpen),sort(openByT{1})));
    capChanges=sum(abs(defaultCaps-capsByT{1}));
    for t=2:numel(openByT)
        lineChanges=lineChanges+numel(setxor(sort(openByT{t-1}),sort(openByT{t})));
        capChanges=capChanges+sum(abs(capsByT{t}-capsByT{t-1}));
    end
end

function writeScheduleCSV(result,sys,schedDef,outDir)
    T=numel(result.openByT);
    TimeStep=(1:T).';
    Mode=strings(T,1);
    OpenLines=strings(T,1);
    CapSteps=strings(T,1);
    Vsp=strings(T,1);
    Vband=strings(T,1);
    InverterAccepted=false(T,1);
    MinV=zeros(T,1); MaxV=zeros(T,1);

    for t=1:T
        if isInWindows(t,schedDef.inverter), Mode(t)='INVERTER'; else, Mode(t)='NETWORK'; end
        OpenLines(t)=string(mat2str(result.openByT{t}));
        CapSteps(t)=string(mat2str(result.capsByT{t}));
        if result.evaluated(t)
            Vsp(t)=string(mat2str(result.VspByT{t},5));
            Vband(t)=string(mat2str(result.VbandByT{t},5));
            InverterAccepted(t)=result.accepted(t);
        else
            Vsp(t)=''; Vband(t)='';
        end
        MinV(t)=min(result.V(t,:));
        MaxV(t)=max(result.V(t,:));
    end

    Tout=table(TimeStep,Mode,OpenLines,CapSteps,Vsp,Vband,InverterAccepted,MinV,MaxV,result.PlossMW, ...
        'VariableNames',{'TimeStep','Mode','OpenLines','CapSteps','Vsp','Vband', ...
        'InverterAccepted','MinV','MaxV','ActiveLossMW'});
    writetable(Tout,fullfile(outDir,sprintf('sensitivity_schedule_%s_%s.csv', ...
        sys.name,schedDef.name)));
end

function T=rowsToTable(rows)
    names={'System','Schedule','Seed','NetworkWindows','InverterWindows', ...
        'JV','TotalViolations','UnderViolations','OverViolations','MinV','MaxV', ...
        'ActiveLoss_kWh','LineStatusChanges','CapStepChanges', ...
        'InverterAccepted','InverterRejected','RuntimeSeconds'};
    T=cell2table(rows,'VariableNames',names);
end

function writeReport(T,cfg,outFile)
    fid=fopen(outFile,'w');
    fprintf(fid,'CONTROL-WINDOW BOUNDARY SENSITIVITY\n');
    fprintf(fid,'===================================\n\n');
    fprintf(fid,'Equal reduced budget used for ALL schedules:\n');
    fprintf(fid,'  network control : %d wolves x %d iterations\n',cfg.netAgents,cfg.netIter);
    fprintf(fid,'  inverter tuning : %d wolves x %d iterations\n',cfg.invAgents,cfg.invIter);
    fprintf(fid,'  PFmin            : %.2f\n',cfg.PFmin);
    fprintf(fid,'  33-bus seed      : %d (reset for each schedule)\n',cfg.seed33);
    fprintf(fid,'  69-bus seed      : %d (reset for each schedule)\n\n',cfg.seed69);
    fprintf(fid,'Interpretation: this experiment isolates +/-30-min schedule-boundary changes.\n');
    fprintf(fid,'Absolute S0 values are NOT intended to reproduce the headline 50x100/40x50 run.\n\n');

    for i=1:height(T)
        fprintf(fid,'%s | %s\n',T.System{i},T.Schedule{i});
        fprintf(fid,'  Network windows : %s\n',T.NetworkWindows{i});
        fprintf(fid,'  Inverter windows: %s\n',T.InverterWindows{i});
        fprintf(fid,'  JV=%.6f | Viol=%d (under=%d over=%d) | V=[%.5f, %.5f]\n', ...
            T.JV(i),T.TotalViolations(i),T.UnderViolations(i),T.OverViolations(i),T.MinV(i),T.MaxV(i));
        fprintf(fid,'  Active loss=%.2f kWh | Line changes=%d | Cap-step changes=%d\n', ...
            T.ActiveLoss_kWh(i),T.LineStatusChanges(i),T.CapStepChanges(i));
        fprintf(fid,'  Inverter accepted/rejected=%d/%d | Runtime=%.1f min\n\n', ...
            T.InverterAccepted(i),T.InverterRejected(i),T.RuntimeSeconds(i)/60);
    end
    fclose(fid);
end

function tf=isInWindows(t,windows)
    tf=false;
    for k=1:size(windows,1)
        if t>=windows(k,1) && t<=windows(k,2)
            tf=true; return;
        end
    end
end

function s=windowsToString(W)
    if size(W,1)==3
        s=sprintf('[%d-%d;%d-%d;%d-%d]',W(1,1),W(1,2),W(2,1),W(2,2),W(3,1),W(3,2));
    elseif size(W,1)==2
        s=sprintf('[%d-%d;%d-%d]',W(1,1),W(1,2),W(2,1),W(2,2));
    else
        parts=strings(size(W,1),1);
        for k=1:size(W,1), parts(k)=sprintf('%d-%d',W(k,1),W(k,2)); end
        s=['[' strjoin(cellstr(parts),';') ']'];
    end
end

function out=ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end


%% ============================= INPUTS ===================================
function systems=makeSystemConfigs(rootDir,cfg)
    %#ok<INUSD>
    systems(1).name='33bus';
    systems(1).caseName='case33bw';
    systems(1).fallbackCaseName='case33bw';
    systems(1).defaultOpen=[33 34 35 36 37];
    systems(1).pvBuses=[21 24 32 17 18];
    systems(1).capBuses=[18 25 33];
    systems(1).pvRatedMW=0.20;
    systems(1).loadFile=locateUniquePreferred(rootDir,'myloadspvs.mat','');
    systems(1).pvFile=systems(1).loadFile;

    systems(2).name='69bus';
    systems(2).caseName='case69modified';
    systems(2).fallbackCaseName='case69';
    systems(2).defaultOpen=[69 70 71 72 73];
    systems(2).pvBuses=[15 18 27 35 46 60];
    systems(2).capBuses=[7 23 30 39 55];
    systems(2).pvRatedMW=1.00;
    systems(2).loadFile=locateUniquePreferred(rootDir,'myloads.mat','');
    systems(2).pvFile=locateUniquePreferred(rootDir,'mypvs.mat','');
end

function [mpc0,Lprof,PVprof,defaultOpen]=loadSystemAndProfiles(sys,cfg)
    try
        mpc0=loadcase(sys.caseName);
    catch
        mpc0=loadcase(sys.fallbackCaseName);
    end
    nBr=size(mpc0.branch,1);
    defaultOpen=sys.defaultOpen(sys.defaultOpen<=nBr);
    if numel(defaultOpen)~=nBr-(size(mpc0.bus,1)-1)
        defaultOpen=find(mpc0.branch(:,11)==0).';
        if isempty(defaultOpen)
            defaultOpen=(nBr-(nBr-(size(mpc0.bus,1)-1))+1):nBr;
        end
    end
    mpc0.branch(:,11)=1;

    if strcmpi(sys.name,'33bus')
        S=load(sys.loadFile);
        if ~isfield(S,'myloads'), error('myloadspvs.mat must contain myloads.'); end
        Lprof=orientProfile(S.myloads,size(mpc0.bus,1),cfg.T);
        if isfield(S,'myPVs'), PVprof=S.myPVs;
        elseif isfield(S,'PV'), PVprof=S.PV;
        else, error('myloadspvs.mat must contain myPVs or PV.'); end
        PVprof=orientProfile(PVprof,numel(sys.pvBuses),cfg.T);
        mx=max(PVprof(:));
        if mx<=0, error('33-bus PV profile has nonpositive maximum.'); end
        PVprof=PVprof./mx.*sys.pvRatedMW;
    else
        L=load(sys.loadFile);
        P=load(sys.pvFile);
        if ~isfield(L,'myloads'), error('myloads.mat must contain myloads.'); end
        Lprof=orientProfile(L.myloads,size(mpc0.bus,1),cfg.T);
        if isfield(P,'PV'), PVprof=P.PV;
        elseif isfield(P,'mypvs'), PVprof=P.mypvs;
        elseif isfield(P,'myPVs'), PVprof=P.myPVs;
        else, error('mypvs.mat must contain PV/mypvs/myPVs.'); end
        PVprof=orientProfile(PVprof,numel(sys.pvBuses),cfg.T);
    end
end

function X=orientProfile(X,nRows,nCols)
    X=double(X);
    if isequal(size(X),[nRows nCols])
        return;
    elseif isequal(size(X),[nCols nRows])
        X=X.';
    else
        error('Profile size %s cannot be oriented to %dx%d.',mat2str(size(X)),nRows,nCols);
    end
end

function root=findProjectRoot(startDir)
    d=startDir;
    while true
        if exist(fullfile(d,'myloadspvs.mat'),'file') && ...
           exist(fullfile(d,'myloads.mat'),'file') && exist(fullfile(d,'mypvs.mat'),'file')
            root=d; return;
        end
        parent=fileparts(d);
        if strcmp(parent,d)
            root=startDir; return;
        end
        d=parent;
    end
end

function f=locateUniquePreferred(rootDir,fileName,preferredToken)
    direct=fullfile(rootDir,fileName);
    if exist(direct,'file')
        f=direct; fprintf('Input resolved: %s\n',f); return;
    end
    D=dir(fullfile(rootDir,'**',fileName));
    D=D(~[D.isdir]);
    if isempty(D), error('Could not find %s under %s.',fileName,rootDir); end
    paths=arrayfun(@(x) fullfile(x.folder,x.name),D,'UniformOutput',false);
    paths=unique(paths);
    if ~isempty(preferredToken)
        hit=contains(paths,preferredToken,'IgnoreCase',true);
        if any(hit), paths=paths(hit); end
    end
    if numel(paths)>1
        % Input matrices should normally have one verified root copy. Prefer the
        % shortest path (usually the project-root file) instead of a stale nested copy.
        [~,idx]=min(cellfun(@length,paths));
        f=paths{idx};
    else
        f=paths{1};
    end
    fprintf('Input resolved: %s\n',f);
end
