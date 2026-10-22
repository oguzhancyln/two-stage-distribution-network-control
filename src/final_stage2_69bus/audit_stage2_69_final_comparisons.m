function audit_stage2_69_final_comparisons()
% Read-only post-optimization IEEE 69-bus comparison, no GWO calls.
% Required: the saved final Stage 2 CSV in results_stage2_reopt_69_1MW,
% MATPOWER case69modified.m + myloads.mat + mypvs.mat in
% ALEXIS_v4_short_names_bundle. Run from project root:
%   octave --no-gui --quiet --eval "audit_stage2_69_final_comparisons()"
%
% Replays the same 96-step discrete schedule with:
%   A: optimized inverter controls, constant-Q capacitors (reference)
%   B: same controls, voltage-dependent capacitor susceptance (BS)
%   C: fixed Vsp=1, Vband=0.03, constant-Q capacitor, acceptance gate
%   D: optimized controls but without acceptance gate
%   E: optimized controls with rating-based Q-limit policy
%
% Keeps all GWO decisions fixed and never overwrites input files.
% Values must be verified by the author against manuscript/Table 15/17.

root=fileparts(mfilename('fullpath'));
bundle=fullfile(root,'ALEXIS_v4_short_names_bundle');
outdir=fullfile(root,'results_stage2_reopt_69_1MW');
addpath(bundle);
if exist('runpf','file')~=2
    mpdir=getenv('MATPOWER_DIR');
    if ~isempty(mpdir) && exist(mpdir,'dir')==7
        addpath(genpath(mpdir));
    end
end
assert(exist('runpf','file')==2 && exist('loadcase','file')==2, 'MATPOWER missing from path.');
filename=fullfile(outdir,'control_schedule_69bus_STAGE2_REOPT.csv');
assert(exist(filename,'file')==2, 'New Stage2 control schedule not found.');
schedule=read_schedule_69(filename);
L=load(fullfile(bundle,'myloads.mat'));
P=load(fullfile(bundle,'mypvs.mat'));
assert(isfield(L,'myloads') && isfield(P,'PV'),'Unexpected MAT field names');
loadProfile=L.myloads; pv=P.PV;
assert(isequal(size(loadProfile),[69,96]) && isequal(size(pv),[6,96]),'Unexpected input profile sizes');
mpc0=loadcase('case69modified');
assert(size(mpc0.bus,1)==69 && size(mpc0.branch,1)==73,'Unexpected MATPOWER feeder');
opt=mpoption('verbose',0,'out.all',0);
cfg.pvBuses=[15 18 27 35 46 60]; cfg.capBuses=[7 23 30 39 55];
cfg.pvRated=1.00; cfg.pfmin=0.90; cfg.capMVAr=0.06;
cfg.vmin=0.95; cfg.vmax=1.05; cfg.lambda=1e5;
modes={'Optimized_ConstantQ','Optimized_ShuntBS','FixedVVC_ConstantQ', ...
       'Optimized_NoGate','RatingCapability_ConstantQ'};
metrics=zeros(numel(modes),10);
byStep=zeros(96,numel(modes)*2+1);
byStep(:,1)=(1:96).';
for c=1:numel(modes)
    vm=zeros(96,69); pLoss=0; qLoss=0; nAccept=0; nReject=0;
    for t=1:96
        s=schedule(t);
        m=mpc0;
        m.branch(:,11)=1;
        m.branch(s.open,11)=0;
        m.bus(:,3)=mpc0.bus(:,3).*loadProfile(:,t);
        m.bus(:,4)=mpc0.bus(:,4).*loadProfile(:,t);
        for k=1:6
            b=cfg.pvBuses(k);
            m.bus(b,3)=m.bus(b,3)-pv(k,t);
        end
        for k=1:5
            b=cfg.capBuses(k);
            if c==2  % fixed susceptance Q(V)=Q_nom*V^2
                m.bus(b,6)=m.bus(b,6)+s.caps(k)*cfg.capMVAr;
            else      % Q constant independently of V
                m.bus(b,4)=m.bus(b,4)-s.caps(k)*cfg.capMVAr;
            end
        end
        r0=runpf(m,opt);
        assert(r0.success==1,sprintf('Pre-action PF failed case=%s step=%d',modes{c},t));
        final=r0;
        if strcmpi(s.mode,'DROOP')
            vpre=r0.bus(:,8);
            jpre=objective(vpre,cfg);
            if c==3
                vsp=ones(1,6); vband=0.03*ones(1,6);
            else
                vsp=s.vsp; vband=s.vband;
            end
            assert(numel(vsp)==6 && numel(vband)==6 && all(vband>0), ...
                   sprintf('Incomplete inverter settings at t=%d',t));
            qmax=min(sqrt(max((cfg.pvRated/cfg.pfmin)^2-pv(:,t).^2,0)), ...
                     pv(:,t)*tan(acos(cfg.pfmin)));
            if c==5
                qmax=min(sqrt(max((cfg.pvRated/cfg.pfmin)^2-pv(:,t).^2,0)), ...
                         cfg.pvRated*tan(acos(cfg.pfmin))*ones(6,1));
            end
            m1=m;
            for k=1:6
                delta=vpre(cfg.pvBuses(k))-vsp(k);
                q=max(-qmax(k),min(qmax(k),-qmax(k)*delta/vband(k)));
                m1.bus(cfg.pvBuses(k),4)=m1.bus(cfg.pvBuses(k),4)-q;
            end
            r1=runpf(m1,opt);
            gate=(c~=4);
            if r1.success && (~gate || objective(r1.bus(:,8),cfg)<=jpre)
                final=r1;
                nAccept=nAccept+1;
            else
                nReject=nReject+1;
            end
        end
        vm(t,:)=final.bus(:,8).';
        lossP=1000*0.25*sum(final.branch(:,14)+final.branch(:,16));
        lossQ=1000*0.25*sum(final.branch(:,15)+final.branch(:,17));
        pLoss=pLoss+lossP; qLoss=qLoss+lossQ;
        byStep(t,2*c)=lossP; byStep(t,2*c+1)=lossQ;
    end
    jv=sum((vm(:)-1).^2);
    violations=nnz(vm(:)<cfg.vmin | vm(:)>cfg.vmax);
    metrics(c,:)=[jv violations min(vm(:)) max(vm(:)) mean(vm(:)) ...
                  pLoss qLoss nAccept nReject nnz(vm(:)<cfg.vmin)];
    fprintf('%-30s JV=%.9f violations=%d Vmin=%.9f Vmax=%.9f P_loss=%.3f kWh Q_loss=%.3f kVArh acc/rej=%d/%d\n', ...
        modes{c},jv,violations,min(vm(:)),max(vm(:)),pLoss,qLoss,nAccept,nReject);
    if c==1
        assert(abs(jv-1.54885382313137)<1e-5 && violations==51 && abs(pLoss-1797.509191829)<0.05, ...
          'REFERENCE REPLAY MISMATCH: stop, investigate inputs and do not report comparisons.');
    end
end
out=fullfile(outdir,'audit_69_new_stage2_comparisons.csv');
fid=fopen(out,'w');assert(fid>0,'Could not create output');
fprintf(fid,'Case,JV,Violations,Vmin,Vmax,Vmean,ActiveLoss_kWh,ReactiveLoss_kVArh,Accepted,Rejected,UnderVoltageCount\n');
for c=1:numel(modes)
    fprintf(fid,'%s,%.12g,%d,%.12g,%.12g,%.12g,%.9f,%.9f,%d,%d,%d\n',modes{c},metrics(c,:));
end
fclose(fid);
out2=fullfile(outdir,'audit_69_new_stage2_comparisons_by_step.csv');
fid=fopen(out2,'w'); assert(fid>0,'Could not create by-step output');
fprintf(fid,'TimeStep');
for c=1:numel(modes)
    fprintf(fid,',%s_Ploss_kWh,%s_Qloss_kVArh',modes{c},modes{c});
end
fprintf(fid,'\n');
for t=1:96
    fprintf(fid,'%d',t);
    for j=2:size(byStep,2), fprintf(fid,',%.12g',byStep(t,j)); end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('CHECKS COMPLETED. Saved:\n%s\n%s\n',out,out2);
end

function rows=read_schedule_69(fn)
 fid=fopen(fn,'r');assert(fid>0,'Could not read schedule');
 header=fgetl(fid); assert(ischar(header)&&~isempty(strfind(header,'TimeStep')),'Bad header');
 rows=repmat(struct('mode','','open',[],'caps',[],'vsp',[],'vband',[]),96,1);
 for t=1:96
   line=fgetl(fid);assert(ischar(line),'Incomplete schedule');
   f=strsplit(line,','); assert(numel(f)==7,'Malformed schedule CSV');
   assert(str2double(f{1})==t,'Time step mismatch');
   rows(t).mode=strtrim(f{2});
   rows(t).open=parse_vec(f{3}); rows(t).caps=parse_vec(f{4});
   rows(t).vsp=parse_vec(f{5});rows(t).vband=parse_vec(f{6});
   assert(numel(rows(t).open)==5 && numel(rows(t).caps)==5,'Missing topology/cap values');
 end
 fclose(fid);
end
function v=parse_vec(s)
 s=strrep(strrep(strtrim(s),'[',''),']','');
 if isempty(s),v=[];else,v=sscanf(s,'%f').';end
end
function J=objective(V,cfg)
 J=sum((V-1).^2)+cfg.lambda*sum(max(cfg.vmin-V,0).^2+max(V-cfg.vmax,0).^2);
end
