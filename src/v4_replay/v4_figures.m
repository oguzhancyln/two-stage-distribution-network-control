%% test01c_consistent_qmax_voltage_envelopes.m
% Rebuilds 96-step voltage-envelope and residual-violation figures from
% the v4 consistent-Qmax safety-gated ablation outputs. Use this AFTER running:
%   run_numeric_ablation_33_69_from_schedules_v4_consistent_qmax.m
%
% Required input files:
%   results_ablation_33_69_v4_consistent_qmax/ablation_results_33bus.mat
%   results_ablation_33_69_v4_consistent_qmax/ablation_results_69bus.mat
%
% Outputs:
%   results_test01c_consistent_qmax_envelopes/
%      test01c_consistent_qmax_summary.csv
%      test01c_consistent_qmax_summary_latex_table.tex
%      fig_33bus_safety_gated_voltage_envelope.png/.fig
%      fig_69bus_safety_gated_voltage_envelope.png/.fig
%      fig_33bus_safety_gated_violation_count_per_step.png/.fig
%      fig_69bus_safety_gated_violation_count_per_step.png/.fig
%
% This script should be used for the final paper figures if the paper uses
% the safety-gated hybrid results from the ablation run.

clear; clc; close all;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end

cfg.Vmin = 0.95;
cfg.Vmax = 1.05;
cfg.inputDir = fullfile(scriptDir, 'results_ablation_33_69_v4_consistent_qmax');
cfg.outputDir = fullfile(scriptDir, 'results_test01c_consistent_qmax_envelopes');
if ~exist(cfg.outputDir, 'dir'), mkdir(cfg.outputDir); end

systems = {'33bus','69bus'};
allRows = table();

for s = 1:numel(systems)
    sysName = systems{s};
    inFile = fullfile(cfg.inputDir, sprintf('ablation_results_%s.mat', sysName));
    if ~exist(inFile, 'file')
        error('Missing %s. Run run_numeric_ablation_33_69_from_schedules_v4_consistent_qmax.m first.', inFile);
    end

    data = load(inFile);
    result = data.result;

    Vno = result.Vcase.NoControl;
    Vhy = result.Vcase.ProposedHybrid;

    % Metrics and summary rows
    mNo = computeMetricsLocal(Vno, cfg.Vmin, cfg.Vmax);
    mHy = computeMetricsLocal(Vhy, cfg.Vmin, cfg.Vmax);

    Tsys = table();
    Tsys.System = {sysName; sysName};
    Tsys.Case = {'No control'; 'Proposed hybrid'};
    Tsys.MinV = [mNo.MinV; mHy.MinV];
    Tsys.MaxV = [mNo.MaxV; mHy.MaxV];
    Tsys.MeanV = [mNo.MeanV; mHy.MeanV];
    Tsys.TotalVoltageDeviation = [mNo.TotalVoltageDeviation; mHy.TotalVoltageDeviation];
    Tsys.UnderVoltageCount = [mNo.UnderVoltageCount; mHy.UnderVoltageCount];
    Tsys.OverVoltageCount = [mNo.OverVoltageCount; mHy.OverVoltageCount];
    Tsys.TotalViolationCount = [mNo.TotalViolationCount; mHy.TotalViolationCount];
    Tsys.DevReductionPercent = [0; 100*(mNo.TotalVoltageDeviation-mHy.TotalVoltageDeviation)/mNo.TotalVoltageDeviation];
    Tsys.ViolationReductionPercent = [0; 100*(mNo.TotalViolationCount-mHy.TotalViolationCount)/max(mNo.TotalViolationCount,eps)];
    allRows = [allRows; Tsys]; %#ok<AGROW>

    fprintf('\n%s safety-gated envelope summary:\n', sysName);
    disp(Tsys);

    makeEnvelopeFigure(sysName, Vno, Vhy, cfg);
    makeViolationFigure(sysName, Vno, Vhy, cfg);
end

writetable(allRows, fullfile(cfg.outputDir, 'test01c_consistent_qmax_summary.csv'));
writeLatexSummaryLocal(allRows, fullfile(cfg.outputDir, 'test01c_consistent_qmax_summary_latex_table.tex'));

fprintf('\nDone. Safety-gated envelope outputs saved under: %s\n', cfg.outputDir);

%% ------------------------- Local functions ----------------------------
function m = computeMetricsLocal(V, Vmin, Vmax)
    under = max(Vmin - V, 0);
    over  = max(V - Vmax, 0);
    m.MinV = min(V(:));
    m.MaxV = max(V(:));
    m.MeanV = mean(V(:));
    m.TotalVoltageDeviation = sum((V(:)-1).^2);
    m.UnderVoltageCount = nnz(under > 0);
    m.OverVoltageCount = nnz(over > 0);
    m.TotalViolationCount = m.UnderVoltageCount + m.OverVoltageCount;
end

function makeEnvelopeFigure(sysName, Vno, Vhy, cfg)
    t = 1:size(Vno,1);
    noMin = min(Vno,[],2); noMax = max(Vno,[],2);
    hyMin = min(Vhy,[],2); hyMax = max(Vhy,[],2);

    fig = figure('Color','w','Position',[100 100 900 600]);
    hold on; grid on; box on;
    plot(t, noMin, '--', 'LineWidth', 1.8);
    plot(t, noMax, '--', 'LineWidth', 1.8);
    plot(t, hyMin, '-', 'LineWidth', 2.2);
    plot(t, hyMax, '-', 'LineWidth', 2.2);
    yline(cfg.Vmin, ':', 'V_{min}', 'LineWidth', 1.2, 'LabelHorizontalAlignment','right');
    yline(cfg.Vmax, ':', 'V_{max}', 'LineWidth', 1.2, 'LabelHorizontalAlignment','right');
    xlabel('Time step');
    ylabel('Voltage magnitude (p.u.)');
    title(sprintf('%s 96-step safety-gated voltage envelope', sysName));
    legend({'No-control min','No-control max','Proposed min','Proposed max'}, 'Location','best');
    ylim([min([noMin;hyMin;cfg.Vmin])-0.02, max([noMax;hyMax;cfg.Vmax])+0.02]);

    base = fullfile(cfg.outputDir, sprintf('fig_%s_safety_gated_voltage_envelope', sysName));
    saveas(fig, [base '.fig']);
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
end

function makeViolationFigure(sysName, Vno, Vhy, cfg)
    t = 1:size(Vno,1);
    noViol = sum(Vno < cfg.Vmin | Vno > cfg.Vmax, 2);
    hyViol = sum(Vhy < cfg.Vmin | Vhy > cfg.Vmax, 2);

    fig = figure('Color','w','Position',[100 100 900 600]);
    hold on; grid on; box on;
    plot(t, noViol, '--', 'LineWidth', 1.8);
    plot(t, hyViol, '-', 'LineWidth', 2.2);
    xlabel('Time step');
    ylabel('Number of violating buses');
    title(sprintf('%s safety-gated residual voltage-limit violations', sysName));
    legend({'No control','Proposed hybrid'}, 'Location','best');
    ylim([0, max([noViol; hyViol])+2]);

    base = fullfile(cfg.outputDir, sprintf('fig_%s_safety_gated_violation_count_per_step', sysName));
    saveas(fig, [base '.fig']);
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
end

function writeLatexSummaryLocal(T, outFile)
    fid = fopen(outFile, 'w');
    if fid < 0, error('Cannot write %s', outFile); end
    fprintf(fid, '\\begin{table}[!t]\n');
    fprintf(fid, '\\caption{Safety-gated voltage-envelope summary for the proposed hybrid strategy.}\n');
    fprintf(fid, '\\label{tab:safety_gated_envelope_summary}\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{llccccc}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'System & Case & Min. V & Max. V & Mean V & $J_V$ & Viol. \\\\ \n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %s & %.4f & %.4f & %.4f & %.4f & %d \\\\ \n', ...
            T.System{i}, T.Case{i}, T.MinV(i), T.MaxV(i), T.MeanV(i), ...
            T.TotalVoltageDeviation(i), T.TotalViolationCount(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
