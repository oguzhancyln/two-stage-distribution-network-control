function audit_stage2_69_losses_octave()
% READ-ONLY 96-step MATPOWER replay for IEEE-69-bus reoptimized Stage 2.
% Run from ALEXIS_PAPER root after placing this file beside stage2_reopt_69_octave.m.
% It NEVER invokes GWO and NEVER alters original schedules or published outputs.
% Requires: bundle with case69modified.m, myloads.mat, mypvs.mat,
% control_schedule_69bus.csv and results_stage2_reopt_69_1MW/*.csv.
% Replays reference, DNR+CAP, legacy 1MW and reoptimized Stage 2.
% Measures losses exactly as stage4_losses_switching.m: (PF+PT), (QF+QT),
% with 15-minute integration and MATPOWER branch columns 14-17.

baseDir = fileparts(mfilename('fullpath'));
bundleDir = fullfile(baseDir, 'ALEXIS_v4_short_names_bundle');
outputDir = fullfile(baseDir, 'results_stage2_reopt_69_1MW');
addpath(bundleDir);
if exist('runpf','file') ~= 2 || exist('loadcase','file') ~= 2
    mpdir = getenv('MATPOWER_DIR');
    if ~isempty(mpdir) && exist(mpdir,'dir') == 7
        addpath(genpath(mpdir));
    end
end
assert(exist('runpf','file') == 2, 'MATPOWER must be on the Octave path.');
oldFile = fullfile(bundleDir, 'control_schedule_69bus.csv');
newFile = fullfile(outputDir, 'control_schedule_69bus_STAGE2_REOPT.csv');
assert(exist(oldFile,'file') == 2 && exist(newFile,'file') == 2, 'Schedules missing.');
old = ar_read_schedule(oldFile);
new = ar_read_schedule(newFile);
for t=1:96
    assert(isequal(old(t).open, new(t).open), 'Open branches changed at time %d.',t);
    assert(isequal(old(t).caps, new(t).caps), 'Capacitor steps changed at time %d.',t);
end
D = load(fullfile(bundleDir,'myloads.mat'));
P = load(fullfile(bundleDir,'mypvs.mat'));
assert(isfield(D,'myloads') && isfield(P,'PV'), 'Unexpected input variable names.');
L = D.myloads; pv = P.PV;
assert(isequal(size(L),[69 96]) && isequal(size(pv),[6 96]), 'Wrong input dimensions.');
mpc0 = loadcase('case69modified');
assert(size(mpc0.bus,1)==69 && size(mpc0.branch,1)==73, 'Unexpected feeder.');
mpopt = mpoption('verbose',0,'out.all',0);
cfg.buses=[15 18 27 35 46 60]; cfg.capBuses=[7 23 30 39 55];
cfg.pf=.90; cfg.pvRatedMW=1; cfg.capStepMVAr=.06;
cfg.vmin=.95; cfg.vmax=1.05; cfg.penalty=1e5; cfg.dth=.25;
defaultOpen=find(mpc0.branch(:,11)==0).';
if numel(defaultOpen)~=5, defaultOpen=69:73; end
names={'Reference','DNR_CAP_only','Legacy_1MW_replay','Reoptimized_Stage2_1MW'};
records=zeros(4,17);
perStep=zeros(96,7); % t, old P, new P, old Q, new Q, old accepted, new accepted
fprintf('Starting READ-ONLY 96-step AC replay. No optimization.\n');
for c=1:4
    V=zeros(96,69); ploss=0; qloss=0;
    accepted=0; rejected=0; nQ=0; sat=0; minPF=inf; minMargin=inf;
    if c==4
        sel=new;
    else
        sel=old;
    end
    openByT=cell(96,1); capsByT=cell(96,1);
    for t=1:96
        if c==1
            op=defaultOpen; caps=zeros(1,5); applyDroop=false;
        else
            op=sel(t).open; caps=sel(t).caps;
            applyDroop=(c>=3 && strcmp(sel(t).mode,'DROOP'));
        end
        openByT{t}=op; capsByT{t}=caps;
        m=mpc0;
        m.branch(:,11)=1; m.branch(op,11)=0;
        m.bus(:,3)=mpc0.bus(:,3).*L(:,t);
        m.bus(:,4)=mpc0.bus(:,4).*L(:,t);
        for k=1:6
            b=cfg.buses(k);
            m.bus(b,3)=m.bus(b,3)-pv(k,t);
        end
        for k=1:5
            m.bus(cfg.capBuses(k),4)=m.bus(cfg.capBuses(k),4)-caps(k)*cfg.capStepMVAr;
        end
        try
            r0=runpf(m,mpopt);
            assert(r0.success==1, 'Pre-action AC PF failed.');
        catch ex
            error('Replay failed at case %s t=%d: %s',names{c},t,ex.message);
        end
        final=r0; qcmd=zeros(6,1); wasAccepted=0;
        if applyDroop
            vpre=abs(r0.bus(:,8));
            Jpre=ar_penalized(vpre,cfg);
            m1=m;
            qmax=min(sqrt(max((cfg.pvRatedMW/cfg.pf)^2-pv(:,t).^2,0)), ...
                     pv(:,t)*tan(acos(cfg.pf)));
            for k=1:6
                d=vpre(cfg.buses(k))-sel(t).vsp(k);
                qcmd(k)=max(-qmax(k),min(qmax(k),-qmax(k)*d/sel(t).vband(k)));
                m1.bus(cfg.buses(k),4)=m1.bus(cfg.buses(k),4)-qcmd(k);
            end
            try
                r1=runpf(m1,mpopt);
                valid=(r1.success==1);
            catch
                valid=false;
            end
            if valid && ar_penalized(abs(r1.bus(:,8)),cfg)<=Jpre
                final=r1; wasAccepted=1; accepted=accepted+1;
                for k=1:6
                    if abs(qcmd(k))>1e-10, nQ=nQ+1; end
                    if qmax(k)>1e-10 && abs(abs(qcmd(k))-qmax(k))<1e-8
                        sat=sat+1;
                    end
                    if pv(k,t)>1e-10
                        minPF=min(minPF,pv(k,t)/hypot(pv(k,t),qcmd(k)));
                        minMargin=min(minMargin,1000*(sqrt(max((cfg.pvRatedMW/cfg.pf)^2-pv(k,t)^2,0))-abs(qcmd(k))));
                    end
                end
            else
                rejected=rejected+1; qcmd(:)=0;
            end
        end
        assert(final.success==1,'Final PF failed');
        V(t,:)=abs(final.bus(:,8)).';
        pStep=sum(final.branch(:,14)+final.branch(:,16))*1000*cfg.dth;
        qStep=sum(final.branch(:,15)+final.branch(:,17))*1000*cfg.dth;
        ploss=ploss+pStep;
        qloss=qloss+qStep;
        if c==3
            perStep(t,1)=t; perStep(t,2)=pStep; perStep(t,4)=qStep;
            perStep(t,6)=wasAccepted;
        elseif c==4
            perStep(t,3)=pStep; perStep(t,5)=qStep;
            perStep(t,7)=wasAccepted;
        end
    end
    jv=sum((V(:)-1).^2);
    viol=nnz(V(:)<cfg.vmin | V(:)>cfg.vmax);
    [lineInit,lineSub,capInit,capSub]=ar_count_switch(openByT,capsByT,defaultOpen);
    if isinf(minPF),minPF=NaN;end
    if isinf(minMargin),minMargin=NaN;end
    records(c,:)=[jv viol min(V(:)) max(V(:)) mean(V(:)) ...
      ploss qloss lineInit lineSub capInit capSub accepted rejected nQ sat minPF minMargin];
    fprintf('%-24s JV=%.12f viol=%d P=%.4f kWh Q=%.4f kVArh accepted/rejected=%d/%d saturated=%d minPF=%.4f headroom=%.3f kVAr\n', ...
      names{c},jv,viol,ploss,qloss,accepted,rejected,sat,minPF,minMargin);
    if c==3
        assert(abs(jv-1.58073212457191)<1e-5 && viol==51, ...
          'Old 1MW result does NOT match published V4; results not finalized.');
    end
    if c==4
        assert(abs(jv-1.54885382313137)<1e-5 && viol==51, ...
          'New 1MW result does NOT match saved Stage2 matrix; results not finalized.');
    end
end
sFile=fullfile(outputDir,'audit_69_full_replay_summary.csv');
fid=fopen(sFile,'w'); assert(fid>0,'Cannot open output summary.');
fprintf(fid,'Method,JV,Violations,Vmin,Vmax,Vmean,ActiveLoss_kWh,ReactiveLoss_kVArh,InitialLineChanges,SubsequentLineChanges,InitialCapSteps,SubsequentCapSteps,Accepted,Rejected,AppliedQCommands,SaturationEvents,MinRealizedPF,MinApparentPowerMargin_kVAr\n');
for c=1:4
    fprintf(fid,'%s,%.12g,%d,%.12g,%.12g,%.12g,%.9f,%.9f,%d,%d,%d,%d,%d,%d,%d,%d,%.9f,%.9f\n', ...
      names{c},records(c,:));
end
fclose(fid);
tFile=fullfile(outputDir,'audit_69_old_new_loss_by_step.csv');
fid=fopen(tFile,'w'); assert(fid>0,'Cannot open output steps.');
fprintf(fid,'TimeStep,LegacyLoss_kWh,NewLoss_kWh,LegacyReactiveLoss_kVArh,NewReactiveLoss_kVArh,LegacyAccepted,NewAccepted\n');
for t=1:96
    fprintf(fid,'%d,%.12g,%.12g,%.12g,%.12g,%d,%d\n',perStep(t,:));
end
fclose(fid);
fprintf('AUDIT COMPLETE. New files:\n%s\n%s\n',sFile,tFile);
end

function rows=ar_read_schedule(fn)
    fid=fopen(fn,'r'); assert(fid>0,'Cannot open schedule.');
    header=fgetl(fid); assert(~isempty(strfind(header,'TimeStep')),'Bad header.');
    rows=repmat(struct('t',0,'mode','','open',[],'caps',[],'vsp',[],'vband',[]),96,1);
    for i=1:96
        line=fgetl(fid); assert(ischar(line),'CSV too short.');
        f=strsplit(line,','); assert(numel(f)==7,'Malformed schedule row.');
        rows(i).t=str2double(f{1}); rows(i).mode=strtrim(f{2});
        rows(i).open=ar_vec(f{3}); rows(i).caps=ar_vec(f{4});
        rows(i).vsp=ar_vec(f{5}); rows(i).vband=ar_vec(f{6});
        assert(rows(i).t==i,'Time step mismatch.');
    end
    fclose(fid);
end
function v=ar_vec(str)
    str=strrep(strrep(strtrim(str),'[',''),']','');
    if isempty(str),v=[];else,v=sscanf(str,'%f').';end
end
function J=ar_penalized(V,cfg)
    J=sum((V-1).^2)+cfg.penalty*sum(max(cfg.vmin-V,0).^2+max(V-cfg.vmax,0).^2);
end
function [li,ls,ci,cs]=ar_count_switch(op,cap,base)
    li=numel(setxor(op{1},base)); ls=0;
    ci=sum(abs(cap{1}));cs=0;
    for t=2:96
        ls=ls+numel(setxor(op{t-1},op{t}));
        cs=cs+sum(abs(cap{t}-cap{t-1}));
    end
end
