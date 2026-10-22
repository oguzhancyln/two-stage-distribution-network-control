%% test03_gwo_pso_de_dnrcap_benchmark.m
% Stage 3 test for the IEEE Access paper:
% Equal-budget GWO / PSO / DE benchmark for the discrete DNR+CAP layer.
%
% IMPORTANT:
%   This benchmark isolates the core discrete search stage (topology
%   reconfiguration + capacitor switching). It does NOT re-optimize the
%   daytime droop layer. The same fitness function, profiles, low-PV windows,
%   radiality check, random radial repair, population size, and iteration
%   budget are used for all three optimizers.
%
% Requirements in MATLAB path/current folder:
%   - MATPOWER with case33bw and case69modified/case69.
%   - 33-bus profiles: myloadspvs.mat with myloads and myPVs/PV.
%   - 69-bus profiles: myloads.mat and mypvs.mat with PV.
%
% Output folder:
%   results_test03_gwo_pso_de_benchmark/
%
% Recommended DEBUG run:
%   cfg.nRuns = 2; cfg.popSize = 12; cfg.maxIter = 15;
%
% Recommended FINAL paper run:
%   cfg.nRuns = 5; cfg.popSize = 30; cfg.maxIter = 60;
%
% The final paper should describe the benchmark as an equal-budget
% comparison for the DNR+CAP layer, not as a full-hybrid benchmark.

clear; clc; close all;
warning('off','all');

%% ---------------- USER SETTINGS ----------------
cfg.outputDir = 'results_test03_gwo_pso_de_benchmark';
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end

% Equal-budget FINAL settings.
cfg.nRuns   = 5;
cfg.popSize = 30;
cfg.maxIter = 60;

cfg.verboseIter = true;
cfg.baseSeed = 31001;

cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.lambda = 1e5;

cfg.maxCapStep = 5;
cfg.capStepMVAr = 60/1000;  % 60 kVAr per step, in MVAr

% Same three DNR+CAP optimization windows used in the manuscript.
cfg.topoWindows = [1 24; 40 60; 73 96];

% PSO parameters.
cfg.psoWmax = 0.90;
cfg.psoWmin = 0.40;
cfg.psoC1 = 1.70;
cfg.psoC2 = 1.70;

% DE parameters.
cfg.deF  = 0.70;
cfg.deCR = 0.90;

methods = {'GWO','PSO','DE'};
systems = makeSystemConfigs();

allRunRows = {};
allAggRows = {};
allConvRows = {};

for sidx = 1:numel(systems)
    sys = systems(sidx);
    fprintf('\n================ Test03 benchmark: %s ================\n', sys.name);

    [mpc0, Lprof, PVprof, defaultOpen] = loadSystemAndProfiles(sys);
    T = size(Lprof,2);
    nb = size(mpc0.bus,1);

    fprintf('Loaded %s: %d buses, %d branches, default open=%s\n', ...
        sys.caseName, nb, size(mpc0.branch,1), mat2str(defaultOpen));
    fprintf('Profiles: load range=[%.3f, %.3f], PV max=%.4f MW/site\n', ...
        min(Lprof(:)), max(Lprof(:)), max(PVprof(:)));

    % No-control radial reference.
    V_no = simulateFullDay(mpc0, Lprof, PVprof, sys, cfg, ...
        repmat(defaultOpen(:)',T,1), zeros(T,numel(sys.capBuses)));
    M_no = computeMetrics(V_no, cfg);

    for midx = 1:numel(methods)
        method = methods{midx};
        fprintf('\n------------ %s / %s ------------\n', sys.name, method);

        for r = 1:cfg.nRuns
            seed = cfg.baseSeed + 10000*sidx + 100*midx + r;
            rng(seed,'twister');
            fprintf('\n--- %s %s run %d/%d | seed=%d ---\n', ...
                sys.name, method, r, cfg.nRuns, seed);

            timerRun = tic;
            openByT = repmat(defaultOpen(:)', T, 1);
            capsByT = zeros(T, numel(sys.capBuses));
            windowObjective = nan(size(cfg.topoWindows,1),1);

            for w = 1:size(cfg.topoWindows,1)
                tRange = cfg.topoWindows(w,1):cfg.topoWindows(w,2);
                label = sprintf('%s %s run %d W%d t=%d..%d', ...
                    sys.name, method, r, w, tRange(1), tRange(end));

                switch upper(method)
                    case 'GWO'
                        [bestOpen,bestCaps,bestObj,conv] = optimizeGWO( ...
                            mpc0,Lprof,PVprof,sys,cfg,tRange,label);
                    case 'PSO'
                        [bestOpen,bestCaps,bestObj,conv] = optimizePSO( ...
                            mpc0,Lprof,PVprof,sys,cfg,tRange,label);
                    case 'DE'
                        [bestOpen,bestCaps,bestObj,conv] = optimizeDE( ...
                            mpc0,Lprof,PVprof,sys,cfg,tRange,label);
                    otherwise
                        error('Unknown method: %s', method);
                end

                windowObjective(w) = bestObj;
                openByT(tRange,:) = repmat(bestOpen(:)',numel(tRange),1);
                capsByT(tRange,:) = repmat(bestCaps(:)',numel(tRange),1);

                % Retain the latest DNR+CAP operating point between windows.
                if w < size(cfg.topoWindows,1)
                    gapRange = (cfg.topoWindows(w,2)+1):(cfg.topoWindows(w+1,1)-1);
                    if ~isempty(gapRange)
                        openByT(gapRange,:) = repmat(bestOpen(:)',numel(gapRange),1);
                        capsByT(gapRange,:) = repmat(bestCaps(:)',numel(gapRange),1);
                    end
                end

                for it = 1:numel(conv)
                    allConvRows(end+1,:) = {sys.name,method,r,w,it,conv(it)}; %#ok<SAGROW>
                end
            end

            % Evaluate the complete 96-step DNR+CAP-only schedule.
            V_sched = simulateFullDay(mpc0,Lprof,PVprof,sys,cfg,openByT,capsByT);
            M = computeMetrics(V_sched,cfg);
            runtime_s = toc(timerRun);

            fprintf('%s %s run %d DONE | J=%.6g | viol=%d | minV=%.5f | maxV=%.5f | runtime=%.1fs\n', ...
                sys.name,method,r,M.TotalVoltageDeviation,M.TotalViolationCount, ...
                M.MinV,M.MaxV,runtime_s);

            row = {sys.name, method, r, seed, ...
                M_no.TotalVoltageDeviation, M_no.TotalViolationCount, ...
                M.TotalVoltageDeviation, M.TotalViolationCount, ...
                M.MinV, M.MaxV, M.MeanV, runtime_s, ...
                windowObjective(1), windowObjective(2), windowObjective(3)};
            allRunRows(end+1,:) = row; %#ok<SAGROW>

            save(fullfile(cfg.outputDir,sprintf('test03_%s_%s_run%02d.mat', ...
                sys.name,lower(method),r)), ...
                'sys','cfg','method','seed','openByT','capsByT','V_no','V_sched', ...
                'M_no','M','windowObjective');
        end
    end
end

%% ---------------- WRITE OUTPUT TABLES ----------------
runNames = {'System','Method','Run','Seed','NoControl_JV','NoControl_Viol', ...
    'Schedule_JV','Schedule_Viol','Schedule_MinV','Schedule_MaxV', ...
    'Schedule_MeanV','RuntimeSeconds','Window1_Obj','Window2_Obj','Window3_Obj'};
T_runs = cell2table(allRunRows,'VariableNames',runNames);
writetable(T_runs,fullfile(cfg.outputDir,'test03_benchmark_all_runs.csv'));

aggRows = {};
for sidx = 1:numel(systems)
    for midx = 1:numel(methods)
        sysName = systems(sidx).name;
        method = methods{midx};
        idx = strcmp(T_runs.System,sysName) & strcmp(T_runs.Method,method);
        D = T_runs(idx,:);
        aggRows(end+1,:) = {sysName,method,height(D), ...
            min(D.Schedule_JV),mean(D.Schedule_JV),std(D.Schedule_JV), ...
            min(D.Schedule_Viol),mean(D.Schedule_Viol),std(D.Schedule_Viol), ...
            mean(D.Schedule_MinV),mean(D.Schedule_MaxV), ...
            median(D.RuntimeSeconds),mean(D.RuntimeSeconds),std(D.RuntimeSeconds)}; %#ok<SAGROW>
    end
end

aggNames = {'System','Method','Runs','Best_JV','Mean_JV','Std_JV', ...
    'Best_Viol','Mean_Viol','Std_Viol','Mean_MinV','Mean_MaxV', ...
    'Median_Runtime_s','Mean_Runtime_s','Std_Runtime_s'};
T_agg = cell2table(aggRows,'VariableNames',aggNames);
writetable(T_agg,fullfile(cfg.outputDir,'test03_benchmark_aggregate.csv'));
writeLatexBenchmarkTable(T_agg,fullfile(cfg.outputDir,'test03_benchmark_aggregate_latex_table.tex'));

if ~isempty(allConvRows)
    convNames = {'System','Method','Run','Window','Iteration','BestObjective'};
    T_conv = cell2table(allConvRows,'VariableNames',convNames);
    writetable(T_conv,fullfile(cfg.outputDir,'test03_benchmark_convergence.csv'));
end

fprintf('\nDone. Test03 outputs saved under: %s\n',cfg.outputDir);

%% ============================ LOCAL FUNCTIONS ============================
function systems = makeSystemConfigs()
    systems(1).name = '33bus';
    systems(1).caseName = 'case33bw';
    systems(1).fallbackCaseName = 'case33bw';
    systems(1).defaultOpen = [33 34 35 36 37];
    systems(1).pvBuses = [21 24 32 17 18];
    systems(1).capBuses = [18 25 33];
    systems(1).pvRatedMW = 0.20;
    systems(1).loadFile = 'myloadspvs.mat';
    systems(1).pvFile = 'myloadspvs.mat';

    systems(2).name = '69bus';
    systems(2).caseName = 'case69modified';
    systems(2).fallbackCaseName = 'case69';
    systems(2).defaultOpen = [69 70 71 72 73];
    systems(2).pvBuses = [15 18 27 35 46 60];
    systems(2).capBuses = [7 23 30 39 55];
    systems(2).pvRatedMW = 1.00;
    systems(2).loadFile = 'myloads.mat';
    systems(2).pvFile = 'mypvs.mat';
end

function [mpc0,Lprof,PVprof,defaultOpen] = loadSystemAndProfiles(sys)
    try
        mpc0 = loadcase(sys.caseName);
    catch
        warning('Could not load %s. Trying %s.',sys.caseName,sys.fallbackCaseName);
        mpc0 = loadcase(sys.fallbackCaseName);
    end

    nBr = size(mpc0.branch,1);
    defaultOpen = sys.defaultOpen(sys.defaultOpen <= nBr);
    if numel(defaultOpen) ~= nBr-(size(mpc0.bus,1)-1)
        defaultOpen = find(mpc0.branch(:,11)==0)';
        if isempty(defaultOpen)
            defaultOpen = (nBr-(nBr-(size(mpc0.bus,1)-1))+1):nBr;
        end
    end
    mpc0.branch(:,11) = 1;

    nb = size(mpc0.bus,1);
    nPV = numel(sys.pvBuses);
    [Lprof,PVprof] = loadProfiles(sys,nb,nPV);

    if strcmpi(sys.name,'33bus')
        mx = max(PVprof(:));
        if mx > 0
            PVprof = PVprof./mx.*sys.pvRatedMW;
        end
    end

    if size(Lprof,2)~=96 || size(PVprof,2)~=96
        error('Profiles must have 96 columns. Found load %s, PV %s.', ...
            mat2str(size(Lprof)),mat2str(size(PVprof)));
    end
    if size(Lprof,1)~=nb
        error('Load profile mismatch for %s: expected %d rows, found %d.', ...
            sys.name,nb,size(Lprof,1));
    end
    if size(PVprof,1)~=nPV
        error('PV profile mismatch for %s: expected %d rows, found %d.', ...
            sys.name,nPV,size(PVprof,1));
    end
end

function [Lprof,PVprof] = loadProfiles(sys,nb,nPV)
    Lprof = []; PVprof = [];
    if strcmpi(sys.name,'33bus') && exist(sys.loadFile,'file')
        S = load(sys.loadFile);
        Lprof = getFieldIfExists(S,{'myloads','loads','L','Lprof'});
        PVprof = getFieldIfExists(S,{'myPVs','PV','mypvs','G','PVprof'});
    elseif strcmpi(sys.name,'69bus')
        if exist(sys.loadFile,'file')
            S = load(sys.loadFile);
            Lprof = getFieldIfExists(S,{'myloads','loads','L','Lprof'});
        end
        if exist(sys.pvFile,'file')
            P = load(sys.pvFile);
            PVprof = getFieldIfExists(P,{'PV','myPVs','mypvs','G','PVprof'});
        end
    end

    if isempty(Lprof) || isempty(PVprof)
        error(['Profiles could not be loaded for %s. This benchmark intentionally ' ...
            'does not generate synthetic profiles.'],sys.name);
    end
    Lprof = orientProfile(Lprof,nb,96);
    PVprof = orientProfile(PVprof,nPV,96);
end

function X = getFieldIfExists(S,names)
    X = [];
    for i = 1:numel(names)
        if isfield(S,names{i}), X = S.(names{i}); return; end
    end
end

function X = orientProfile(X,nRows,nCols)
    X = double(X);
    if isequal(size(X),[nRows nCols])
        return;
    elseif isequal(size(X),[nCols nRows])
        X = X';
    else
        error('Profile size %s cannot be oriented to %dx%d.', ...
            mat2str(size(X)),nRows,nCols);
    end
end

function Vmat = simulateFullDay(mpc0,Lprof,PVprof,sys,cfg,openByT,capsByT)
    T = size(Lprof,2);
    nb = size(mpc0.bus,1);
    Vmat = nan(T,nb);
    for t = 1:T
        Vmat(t,:) = runOneTimeStep(mpc0,Lprof,PVprof,sys,cfg,t, ...
            openByT(t,:),capsByT(t,:));
    end
end

function V = runOneTimeStep(mpc0,Lprof,PVprof,sys,cfg,t,openLines,capSteps)
    mpc = mpc0;
    mpc.branch(:,11) = 1;
    openLines = unique(round(openLines(:)'));
    openLines = openLines(openLines>=1 & openLines<=size(mpc.branch,1));
    mpc.branch(openLines,11) = 0;

    P0 = mpc.bus(:,3); Q0 = mpc.bus(:,4);
    mpc.bus(:,3) = P0.*Lprof(:,t);
    mpc.bus(:,4) = Q0.*Lprof(:,t);

    for k = 1:numel(sys.pvBuses)
        mpc.bus(sys.pvBuses(k),3) = mpc.bus(sys.pvBuses(k),3)-PVprof(k,t);
    end

    capSteps = round(capSteps(:)');
    for c = 1:numel(sys.capBuses)
        mpc.bus(sys.capBuses(c),4) = mpc.bus(sys.capBuses(c),4)- ...
            capSteps(c)*cfg.capStepMVAr;
    end

    mpopt = mpoption('verbose',0,'out.all',0);
    r = runpf(mpc,mpopt);
    if ~r.success
        V = 0.5*ones(1,size(mpc0.bus,1));
    else
        V = abs(r.bus(:,8))';
    end
end

function score = fitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps)
    score = 0;
    for t = tRange
        V = runOneTimeStep(mpc0,Lprof,PVprof,sys,cfg,t,openLines,caps);
        score = score + timeStepObjective(V,cfg);
    end
end

function J = timeStepObjective(V,cfg)
    under = max(cfg.Vmin-V,0);
    over  = max(V-cfg.Vmax,0);
    J = sum((V-1).^2)+cfg.lambda*sum(under.^2+over.^2);
end

function M = computeMetrics(Vmat,cfg)
    under = Vmat<cfg.Vmin;
    over  = Vmat>cfg.Vmax;
    M.MinV = min(Vmat(:));
    M.MaxV = max(Vmat(:));
    M.MeanV = mean(Vmat(:));
    M.TotalVoltageDeviation = sum((Vmat(:)-1).^2);
    M.TotalViolationCount = nnz(under|over);
end

function [lb,ub,nOpen,nCap,dim] = decisionBounds(mpc0,sys,cfg)
    nb = size(mpc0.bus,1);
    nBr = size(mpc0.branch,1);
    nOpen = nBr-(nb-1);
    nCap = numel(sys.capBuses);
    dim = nOpen+nCap;
    lb = [ones(1,nOpen),zeros(1,nCap)];
    ub = [nBr*ones(1,nOpen),cfg.maxCapStep*ones(1,nCap)];
end

function X = initialPopulation(mpc0,sys,cfg)
    [~,~,nOpen,nCap,dim] = decisionBounds(mpc0,sys,cfg);
    X = zeros(cfg.popSize,dim);
    for a = 1:cfg.popSize
        X(a,1:nOpen) = randomRadialOpen(mpc0,nOpen);
        X(a,nOpen+1:end) = randi([0 cfg.maxCapStep],1,nCap);
    end
end

function [openLines,caps,xDecoded] = decode(x,mpc0,sys,cfg)
    [~,~,nOpen,nCap,~] = decisionBounds(mpc0,sys,cfg);
    nBr = size(mpc0.branch,1);
    openLines = unique(round(x(1:nOpen)),'stable');
    openLines = openLines(openLines>=1 & openLines<=nBr);
    if numel(openLines)~=nOpen || ~isRadialAfterOpening(mpc0,openLines)
        openLines = randomRadialOpen(mpc0,nOpen);
    end
    caps = round(x(nOpen+1:nOpen+nCap));
    caps = max(min(caps,cfg.maxCapStep),0);
    xDecoded = [openLines caps];
end

function openLines = randomRadialOpen(mpc0,nOpen)
    nBr = size(mpc0.branch,1);
    for tries = 1:10000
        openLines = sort(randperm(nBr,nOpen));
        if isRadialAfterOpening(mpc0,openLines), return; end
    end
    error('Could not sample a radial opening. Check topology.');
end

function ok = isRadialAfterOpening(mpc0,openLines)
    nb = size(mpc0.bus,1);
    nBr = size(mpc0.branch,1);
    status = ones(nBr,1); status(openLines)=0;
    closed = mpc0.branch(status==1,1:2);

    if size(closed,1)~=nb-1
        ok = false; return;
    end

    visited = false(nb,1);
    q = 1; visited(1)=true;
    while ~isempty(q)
        u = q(1); q(1)=[];
        idx = closed(:,1)==u | closed(:,2)==u;
        edges = closed(idx,:);
        for e = 1:size(edges,1)
            v = edges(e,1)+edges(e,2)-u;
            if v>=1 && v<=nb && ~visited(v)
                visited(v)=true;
                q(end+1)=v; %#ok<AGROW>
            end
        end
    end
    ok = all(visited);
end

function [bestOpen,bestCaps,bestObj,conv] = optimizeGWO(mpc0,Lprof,PVprof,sys,cfg,tRange,label)
    [lb,ub,~,~,~] = decisionBounds(mpc0,sys,cfg);
    X = initialPopulation(mpc0,sys,cfg);

    alpha = X(1,:); beta = X(1,:); delta = X(1,:);
    alphaScore = inf; betaScore = inf; deltaScore = inf;
    conv = nan(cfg.maxIter,1);

    for it = 1:cfg.maxIter
        for a = 1:cfg.popSize
            [openLines,caps,xDec] = decode(X(a,:),mpc0,sys,cfg);
            X(a,:) = xDec;
            score = fitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps);

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

        aGWO = 2-it*(2/cfg.maxIter);
        X = updateGWOPositions(X,alpha,beta,delta,aGWO,lb,ub);
        conv(it)=alphaScore;
        printProgress(cfg,label,it,alphaScore);
    end
    [bestOpen,bestCaps] = decode(alpha,mpc0,sys,cfg);
    bestObj = alphaScore;
end

function Xnew = updateGWOPositions(X,alpha,beta,delta,aGWO,lb,ub)
    Xnew = X;
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            r1=rand; r2=rand; A1=2*aGWO*r1-aGWO; C1=2*r2;
            X1=alpha(j)-A1*abs(C1*alpha(j)-X(i,j));
            r1=rand; r2=rand; A2=2*aGWO*r1-aGWO; C2=2*r2;
            X2=beta(j)-A2*abs(C2*beta(j)-X(i,j));
            r1=rand; r2=rand; A3=2*aGWO*r1-aGWO; C3=2*r2;
            X3=delta(j)-A3*abs(C3*delta(j)-X(i,j));
            Xnew(i,j)=(X1+X2+X3)/3;
        end
    end
    Xnew=max(min(Xnew,ub),lb);
end

function [bestOpen,bestCaps,bestObj,conv] = optimizePSO(mpc0,Lprof,PVprof,sys,cfg,tRange,label)
    [lb,ub,~,~,dim] = decisionBounds(mpc0,sys,cfg);
    X = initialPopulation(mpc0,sys,cfg);
    V = zeros(cfg.popSize,dim);
    pBest = X;
    pBestScore = inf(cfg.popSize,1);
    gBest = X(1,:);
    gBestScore = inf;
    conv = nan(cfg.maxIter,1);

    vmax = 0.20*(ub-lb);

    for it = 1:cfg.maxIter
        for a = 1:cfg.popSize
            [openLines,caps,xDec] = decode(X(a,:),mpc0,sys,cfg);
            X(a,:) = xDec;
            score = fitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps);

            if score<pBestScore(a)
                pBestScore(a)=score;
                pBest(a,:)=X(a,:);
            end
            if score<gBestScore
                gBestScore=score;
                gBest=X(a,:);
            end
        end

        w = cfg.psoWmax-(cfg.psoWmax-cfg.psoWmin)*(it-1)/max(cfg.maxIter-1,1);
        for a = 1:cfg.popSize
            r1 = rand(1,dim); r2 = rand(1,dim);
            V(a,:) = w*V(a,:)+cfg.psoC1*r1.*(pBest(a,:)-X(a,:)) + ...
                cfg.psoC2*r2.*(gBest-X(a,:));
            V(a,:) = max(min(V(a,:),vmax),-vmax);
            X(a,:) = X(a,:)+V(a,:);
        end
        X=max(min(X,ub),lb);
        conv(it)=gBestScore;
        printProgress(cfg,label,it,gBestScore);
    end
    [bestOpen,bestCaps] = decode(gBest,mpc0,sys,cfg);
    bestObj = gBestScore;
end

function [bestOpen,bestCaps,bestObj,conv] = optimizeDE(mpc0,Lprof,PVprof,sys,cfg,tRange,label)
    [lb,ub,~,~,dim] = decisionBounds(mpc0,sys,cfg);
    X = initialPopulation(mpc0,sys,cfg);
    scores = inf(cfg.popSize,1);

    for a = 1:cfg.popSize
        [openLines,caps,xDec] = decode(X(a,:),mpc0,sys,cfg);
        X(a,:) = xDec;
        scores(a) = fitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps);
    end

    [bestObj,bestIdx] = min(scores);
    bestX = X(bestIdx,:);
    conv = nan(cfg.maxIter,1);

    for it = 1:cfg.maxIter
        for i = 1:cfg.popSize
            idxPool = setdiff(1:cfg.popSize,i);
            rr = idxPool(randperm(numel(idxPool),3));
            mutant = X(rr(1),:)+cfg.deF*(X(rr(2),:)-X(rr(3),:));
            mutant=max(min(mutant,ub),lb);

            trial = X(i,:);
            jrand = randi(dim);
            mask = rand(1,dim)<=cfg.deCR;
            mask(jrand)=true;
            trial(mask)=mutant(mask);

            [openLines,caps,trialDec] = decode(trial,mpc0,sys,cfg);
            trialScore = fitness(mpc0,Lprof,PVprof,sys,cfg,tRange,openLines,caps);

            if trialScore<scores(i)
                X(i,:)=trialDec;
                scores(i)=trialScore;
                if trialScore<bestObj
                    bestObj=trialScore;
                    bestX=trialDec;
                end
            end
        end
        conv(it)=bestObj;
        printProgress(cfg,label,it,bestObj);
    end
    [bestOpen,bestCaps] = decode(bestX,mpc0,sys,cfg);
end

function printProgress(cfg,label,it,bestObj)
    if cfg.verboseIter && (it==1 || mod(it,10)==0 || it==cfg.maxIter)
        fprintf('[%s] iter %3d/%3d | best=%12.6g\n',label,it,cfg.maxIter,bestObj);
    end
end

function writeLatexBenchmarkTable(T,filePath)
    fid = fopen(filePath,'w');
    if fid<0, error('Cannot write %s',filePath); end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid,'%% Auto-generated by test03_gwo_pso_de_dnrcap_benchmark.m\n');
    fprintf(fid,'\\begin{table}[!t]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\caption{Equal-budget optimizer benchmark for the DNR+CAP layer.}\n');
    fprintf(fid,'\\label{tab:dnrcap_optimizer_benchmark}\n');
    fprintf(fid,'\\begin{tabular}{llrrrrr}\n');
    fprintf(fid,'\\hline\n');
    fprintf(fid,'System & Method & Best $J_V$ & Mean $J_V$ & Std. $J_V$ & Mean viol. & Median runtime (s) \\\\ \\hline\n');

    for i=1:height(T)
        fprintf(fid,'%s & %s & %.4f & %.4f & %.4f & %.1f & %.1f \\\\ \n', ...
            T.System{i},T.Method{i},T.Best_JV(i),T.Mean_JV(i),T.Std_JV(i), ...
            T.Mean_Viol(i),T.Median_Runtime_s(i));
    end

    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{table}\n');
end
