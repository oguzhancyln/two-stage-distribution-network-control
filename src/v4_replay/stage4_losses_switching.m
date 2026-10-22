%% stage4_losses_switching.m
% Replay-based loss and switching-effort diagnostic for the selected V4
% consistent-Qmax schedules. No optimization is performed.
%
% Run from the same folder as:
%   control_schedule_33bus.csv, control_schedule_69bus.csv
%   myloadspvs.mat, myloads.mat, mypvs.mat
%
% Requirements:
%   MATPOWER on path; case33bw and case69modified available.
%
% Outputs:
%   results_stage4_losses_switching/stage4_loss_switch_summary.csv
%   results_stage4_losses_switching/stage4_loss_switch_summary_latex_table.tex
%
% Definitions:
% - Active/reactive loss energies are obtained from MATPOWER branch flows.
% - Initial line-state changes compare the selected t=1 topology with the
%   conventional default topology.
% - Subsequent line-state changes count status flips between successive
%   time steps. Each physical switch opening or closing counts as one.
% - Initial capacitor-step changes compare t=1 steps with an all-zero
%   reference. Subsequent capacitor-step changes sum absolute step changes.
% - Proposed hybrid uses the same V4 safety gate and conservative Qmax model.

clear; clc;
warning('off','all');

thisFile = mfilename('fullpath');
if isempty(thisFile), scriptDir = pwd; else, scriptDir = fileparts(thisFile); end
addpath(scriptDir);

cfg.T = 96;
cfg.dtHours = 0.25;
cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambdaUnder = 1e5;
cfg.lambdaOver = 1e5;
cfg.PFmin = 0.90;
cfg.outputDir = fullfile(scriptDir,'results_stage4_losses_switching');
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end

case33 = struct('name','33bus','caseName','case33bw', ...
    'scheduleFile',fullfile(scriptDir,'control_schedule_33bus.csv'), ...
    'profileType','combined33','profileFile',fullfile(scriptDir,'myloadspvs.mat'), ...
    'pvBuses',[21 24 32 17 18],'capBuses',[18 25 33], ...
    'capStepMVAr',60/1000,'pvMode','normalized_to_rating','pvRatedMW',0.20);

case69 = struct('name','69bus','caseName','case69modified', ...
    'scheduleFile',fullfile(scriptDir,'control_schedule_69bus.csv'), ...
    'profileType','split69','loadFile',fullfile(scriptDir,'myloads.mat'), ...
    'pvFile',fullfile(scriptDir,'mypvs.mat'), ...
    'pvBuses',[15 18 27 35 46 60],'capBuses',[7 23 30 39 55], ...
    'capStepMVAr',60/1000,'pvMode','MW','pvRatedMW',1.00);

R33 = runSystem(case33,cfg);
R69 = runSystem(case69,cfg);
S = [R33; R69];

writetable(S,fullfile(cfg.outputDir,'stage4_loss_switch_summary.csv'));
writeLatexTable(S,fullfile(cfg.outputDir,'stage4_loss_switch_summary_latex_table.tex'));

disp(S);
fprintf('\nDone. Stage-4 diagnostic outputs saved under: %s\n',cfg.outputDir);

%% ============================= FUNCTIONS =============================
function S = runSystem(sys,cfg)
    fprintf('\n================ %s replay diagnostic ================\n',sys.name);
    mpc0=loadcase(sys.caseName);
    nb=size(mpc0.bus,1); nbr=size(mpc0.branch,1);
    nOpen=nbr-(nb-1);
    baseOpen=determineDefaultOpenLines(mpc0,nOpen);
    [loadMult,Ppv]=loadProfiles(sys,cfg.T,nb,numel(sys.pvBuses));
    sched=readSchedule(sys.scheduleFile,cfg.T);

    cases={'No control','DNR+CAP only','Proposed hybrid'};
    rows=cell(numel(cases),12);

    for c=1:numel(cases)
        name=cases{c};
        fprintf('Evaluating %s / %s ...\n',sys.name,name);
        EpkWh=0; EqkVArh=0; accepted=0; rejected=0;
        openByT=cell(cfg.T,1); capByT=cell(cfg.T,1);

        for t=1:cfg.T
            row=sched(t);
            switch name
                case 'No control'
                    openLines=baseOpen;
                    capSteps=zeros(1,numel(sys.capBuses));
                    useDroop=false; Vsp=[]; Vband=[];
                case 'DNR+CAP only'
                    openLines=row.OpenLines; if isempty(openLines), openLines=baseOpen; end
                    capSteps=row.CapSteps; if isempty(capSteps), capSteps=zeros(1,numel(sys.capBuses)); end
                    useDroop=false; Vsp=[]; Vband=[];
                case 'Proposed hybrid'
                    openLines=row.OpenLines; if isempty(openLines), openLines=baseOpen; end
                    capSteps=row.CapSteps; if isempty(capSteps), capSteps=zeros(1,numel(sys.capBuses)); end
                    useDroop=row.HasDroop;
                    Vsp=row.Vsp; Vband=row.Vband;
            end

            openByT{t}=sort(openLines(:).');
            capByT{t}=capSteps(:).';

            [rFinal,wasAccepted,wasRejected]=runSelectedPF(mpc0,openLines,capSteps,useDroop,Vsp,Vband,t,loadMult,Ppv,sys,cfg);
            accepted=accepted+wasAccepted;
            rejected=rejected+wasRejected;

            if ~rFinal.success
                error('%s / %s / t=%d: power flow did not converge.',sys.name,name,t);
            end

            % MATPOWER result columns: PF=14, QF=15, PT=16, QT=17.
            PlossMW=sum(rFinal.branch(:,14)+rFinal.branch(:,16));
            QlossMVAr=sum(rFinal.branch(:,15)+rFinal.branch(:,17));
            EpkWh=EpkWh+1000*PlossMW*cfg.dtHours;
            EqkVArh=EqkVArh+1000*QlossMVAr*cfg.dtHours;
        end

        [lineInit,lineSub]=switchEffort(openByT,baseOpen);
        [capInit,capSub]=capEffort(capByT,numel(sys.capBuses));

        rows(c,:)={sys.name,name,EpkWh,EqkVArh,lineInit,lineSub,lineInit+lineSub, ...
            capInit,capSub,capInit+capSub,accepted,rejected};
    end

    vars={'System','Case','ActiveLossEnergy_kWh','ReactiveLossEnergy_kVArh', ...
        'InitialLineStateChanges','SubsequentLineStateChanges','TotalLineStateChanges', ...
        'InitialCapStepChanges','SubsequentCapStepChanges','TotalCapStepChanges', ...
        'DroopAcceptedSteps','DroopRejectedSteps'};
    S=cell2table(rows,'VariableNames',vars);
end

function [rFinal,accepted,rejected]=runSelectedPF(mpc0,openLines,capSteps,useDroop,Vsp,Vband,t,loadMult,Ppv,sys,cfg)
    mpopt=mpoption('verbose',0,'out.all',0);
    m0=prepareCase(mpc0,openLines,capSteps,t,loadMult,Ppv,sys);
    r0=runpf(m0,mpopt);
    accepted=0; rejected=0;
    if ~useDroop || ~r0.success
        rFinal=r0; return;
    end

    V0=abs(r0.bus(:,8)).';
    J0=objectiveVoltage(V0,cfg);
    Qmax=computeQmax(Ppv(:,t),mpc0.baseMVA,cfg,sys);
    m1=m0;
    for k=1:numel(sys.pvBuses)
        b=sys.pvBuses(k);
        Qpu=pvControl(V0(b),Vsp(k),Vband(k),Qmax(k));
        m1.bus(b,4)=m1.bus(b,4)-Qpu*mpc0.baseMVA;
    end
    r1=runpf(m1,mpopt);
    if ~r1.success
        rFinal=r0; rejected=1; return;
    end

    V1=abs(r1.bus(:,8)).';
    J1=objectiveVoltage(V1,cfg);
    if J1<=J0
        rFinal=r1; accepted=1;
    else
        rFinal=r0; rejected=1;
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

function qmax=computeQmax(PMW,baseMVA,cfg,sys)
    SmaxMW=sys.pvRatedMW/cfg.PFmin;
    qCircle=sqrt(max(SmaxMW.^2-PMW.^2,0));
    qPf=PMW.*tan(acos(cfg.PFmin));
    qmax=min(qCircle,qPf)./baseMVA;
end

function q=pvControl(V,Vsp,Vband,Qmax)
    dV=V-Vsp;
    q=-(Qmax/Vband)*dV;
    if dV<=-Vband, q=Qmax; elseif dV>=Vband, q=-Qmax; end
end

function J=objectiveVoltage(V,cfg)
    under=max(cfg.Vmin-V,0); over=max(V-cfg.Vmax,0);
    J=sum((V-1).^2)+cfg.lambdaUnder*sum(under.^2)+cfg.lambdaOver*sum(over.^2);
end

function [initCount,subCount]=switchEffort(openByT,baseOpen)
    initCount=numel(setxor(sort(baseOpen(:).'),sort(openByT{1}(:).')));
    subCount=0;
    for t=2:numel(openByT)
        subCount=subCount+numel(setxor(sort(openByT{t-1}(:).'),sort(openByT{t}(:).')));
    end
end

function [initCount,subCount]=capEffort(capByT,nCap)
    z=zeros(1,nCap);
    initCount=sum(abs(capByT{1}-z));
    subCount=0;
    for t=2:numel(capByT)
        subCount=subCount+sum(abs(capByT{t}-capByT{t-1}));
    end
end

function [L,P]=loadProfiles(sys,T,nb,nPV)
    switch sys.profileType
        case 'combined33'
            D=load(sys.profileFile);
            L=D.myloads; P=D.myPVs;
        case 'split69'
            A=load(sys.loadFile); B=load(sys.pvFile);
            L=A.myloads; P=B.PV;
    end
    if size(L,1)~=nb || size(L,2)~=T, error('%s: invalid load matrix size.',sys.name); end
    if size(P,1)~=nPV || size(P,2)~=T, error('%s: invalid PV matrix size.',sys.name); end
    if strcmp(sys.pvMode,'normalized_to_rating')
        mx=max(P(:)); if mx>0, P=P./mx.*sys.pvRatedMW; end
    end
end

function sched=readSchedule(filename,T)
    TBL=readtable(filename,'TextType','string');
    sched=repmat(struct('OpenLines',[],'CapSteps',[],'Vsp',[],'Vband',[],'HasDroop',false),T,1);
    for i=1:height(TBL)
        t=TBL.TimeStep(i); if t<1 || t>T, continue; end
        sched(t).OpenLines=parseVector(TBL.OpenLines(i));
        sched(t).CapSteps=parseVector(TBL.CapSteps(i));
        sched(t).Vsp=parseVector(TBL.Vsp(i));
        sched(t).Vband=parseVector(TBL.Vband(i));
        sched(t).HasDroop=~isempty(sched(t).Vsp) && ~isempty(sched(t).Vband);
    end
end

function v=parseVector(x)
    if ismissing(x), v=[]; return; end
    s=strtrim(char(x));
    if isempty(s) || strcmp(s,'[]'), v=[]; return; end
    s=strrep(s,'[',''); s=strrep(s,']',''); s=strrep(s,',',' ');
    v=sscanf(s,'%f').';
end

function openLines=determineDefaultOpenLines(mpc,nOpen)
    existing=find(mpc.branch(:,11)==0).';
    if numel(existing)==nOpen, openLines=existing;
    else, openLines=(size(mpc.branch,1)-nOpen+1):size(mpc.branch,1);
    end
end

function writeLatexTable(S,filename)
    fid=fopen(filename,'w');
    fprintf(fid,'%% Auto-generated by stage4_losses_switching.m\n');
    fprintf(fid,'\\begin{table*}[!t]\\n\\centering\\n');
    fprintf(fid,'\\caption{Replay-based loss and switching-effort diagnostic for the selected schedules.}\\n');
    fprintf(fid,'\\label{tab:loss_switch_diag}\\n\\scriptsize\\n');
    fprintf(fid,'\\begin{tabular}{llrrrrrr}\\n\\hline\\n');
    fprintf(fid,'System & Case & $E_{P,loss}$ (kWh) & $E_{Q,loss}$ (kVArh) & Line init. & Line subsequent & Cap init. & Cap subsequent \\\\\\\\ \\hline\\n');
    for i=1:height(S)
        fprintf(fid,'%s & %s & %.2f & %.2f & %d & %d & %d & %d \\\\\\\\ \\n', ...
            S.System{i},S.Case{i},S.ActiveLossEnergy_kWh(i),S.ReactiveLossEnergy_kVArh(i), ...
            S.InitialLineStateChanges(i),S.SubsequentLineStateChanges(i), ...
            S.InitialCapStepChanges(i),S.SubsequentCapStepChanges(i));
    end
    fprintf(fid,'\\hline\\n\\end{tabular}\\n\\end{table*}\\n');
    fclose(fid);
end
