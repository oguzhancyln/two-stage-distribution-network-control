%% v4_figures.m
% Rebuilds the 96-step voltage-profile and voltage-limit-violation figures
% from the saved v4 consistent-Qmax results.
%
% IMPORTANT:
% - This script does NOT rerun optimization.
% - Run it only after the consistent-Qmax result files have been generated.
%
% Required input files:
%   results_ablation_33_69_v4_consistent_qmax/ablation_results_33bus.mat
%   results_ablation_33_69_v4_consistent_qmax/ablation_results_69bus.mat
%
% Main outputs:
%   results_test01c_consistent_qmax_envelopes/
%      test01c_consistent_qmax_summary.csv
%      test01c_consistent_qmax_summary_latex_table.tex
%      fig_33bus_daily_voltage_profiles.png/.fig
%      fig_69bus_daily_voltage_profiles.png/.fig
%      fig_33bus_voltage_limit_violations.png/.fig
%      fig_69bus_voltage_limit_violations.png/.fig
%
% Backward-compatible copies are also written using the former filenames,
% so the current LaTeX file does not have to be changed immediately.

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

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

systems = {'33bus', '69bus'};
allRows = table();

for s = 1:numel(systems)
    sysName = systems{s};

    inFile = fullfile(cfg.inputDir, sprintf('ablation_results_%s.mat', sysName));
    if ~exist(inFile, 'file')
        error(['Missing input file:\n  %s\n\n', ...
               'Run the consistent-Qmax analysis first.'], inFile);
    end

    data = load(inFile);

    if ~isfield(data, 'result') || ~isfield(data.result, 'Vcase') || ...
       ~isfield(data.result.Vcase, 'NoControl') || ...
       ~isfield(data.result.Vcase, 'ProposedHybrid')
        error('Unexpected MAT-file structure in:\n  %s', inFile);
    end

    result = data.result;

    % Existing internal field names are retained because they come from the
    % saved result files. Visible labels in figures and tables are updated.
    Vref = result.Vcase.NoControl;
    Vcoord = result.Vcase.ProposedHybrid;

    % Metrics and summary rows
    mRef = computeMetricsLocal(Vref, cfg.Vmin, cfg.Vmax);
    mCoord = computeMetricsLocal(Vcoord, cfg.Vmin, cfg.Vmax);

    Tsys = table();
    Tsys.System = {sysName; sysName};
    Tsys.Case = {'Reference operation'; 'Coordinated strategy'};
    Tsys.MinV = [mRef.MinV; mCoord.MinV];
    Tsys.MaxV = [mRef.MaxV; mCoord.MaxV];
    Tsys.MeanV = [mRef.MeanV; mCoord.MeanV];
    Tsys.TotalVoltageDeviation = ...
        [mRef.TotalVoltageDeviation; mCoord.TotalVoltageDeviation];
    Tsys.UnderVoltageCount = ...
        [mRef.UnderVoltageCount; mCoord.UnderVoltageCount];
    Tsys.OverVoltageCount = ...
        [mRef.OverVoltageCount; mCoord.OverVoltageCount];
    Tsys.TotalViolationCount = ...
        [mRef.TotalViolationCount; mCoord.TotalViolationCount];
    Tsys.DevReductionPercent = ...
        [0; 100 * (mRef.TotalVoltageDeviation - ...
                   mCoord.TotalVoltageDeviation) / ...
                   max(mRef.TotalVoltageDeviation, eps)];
    Tsys.ViolationReductionPercent = ...
        [0; 100 * (mRef.TotalViolationCount - ...
                   mCoord.TotalViolationCount) / ...
                   max(mRef.TotalViolationCount, eps)];

    allRows = [allRows; Tsys]; %#ok<AGROW>

    fprintf('\n%s daily voltage-quality summary:\n', ...
        systemLabelLocal(sysName));
    disp(Tsys);

    makeVoltageProfileFigureLocal(sysName, Vref, Vcoord, cfg);
    makeViolationFigureLocal(sysName, Vref, Vcoord, cfg);
end

writetable(allRows, ...
    fullfile(cfg.outputDir, 'test01c_consistent_qmax_summary.csv'));

writeLatexSummaryLocal(allRows, ...
    fullfile(cfg.outputDir, ...
    'test01c_consistent_qmax_summary_latex_table.tex'));

fprintf('\nDone. Updated figures and summary files were saved under:\n  %s\n', ...
    cfg.outputDir);

%% ------------------------- Local functions ----------------------------

function m = computeMetricsLocal(V, Vmin, Vmax)
    under = max(Vmin - V, 0);
    over  = max(V - Vmax, 0);

    m.MinV = min(V(:));
    m.MaxV = max(V(:));
    m.MeanV = mean(V(:));
    m.TotalVoltageDeviation = sum((V(:) - 1).^2);
    m.UnderVoltageCount = nnz(under > 0);
    m.OverVoltageCount = nnz(over > 0);
    m.TotalViolationCount = ...
        m.UnderVoltageCount + m.OverVoltageCount;
end

function makeVoltageProfileFigureLocal(sysName, Vref, Vcoord, cfg)
    t = 1:size(Vref, 1);

    refMin = min(Vref, [], 2);
    refMax = max(Vref, [], 2);
    coordMin = min(Vcoord, [], 2);
    coordMax = max(Vcoord, [], 2);

    fig = figure( ...
        'Color', 'w', ...
        'Position', [100 100 900 600], ...
        'MenuBar', 'none', ...
        'ToolBar', 'none');

    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    plot(ax, t, refMin, '--', 'LineWidth', 1.8);
    plot(ax, t, refMax, '--', 'LineWidth', 1.8);
    plot(ax, t, coordMin, '-', 'LineWidth', 2.2);
    plot(ax, t, coordMax, '-', 'LineWidth', 2.2);

    yline(ax, cfg.Vmin, ':', 'V_{min}', ...
        'LineWidth', 1.2, ...
        'LabelHorizontalAlignment', 'right');

    yline(ax, cfg.Vmax, ':', 'V_{max}', ...
        'LineWidth', 1.2, ...
        'LabelHorizontalAlignment', 'right');

    xlabel(ax, 'Time step');
    ylabel(ax, 'Voltage magnitude (p.u.)');

    title(ax, sprintf('%s minimum and maximum voltage profiles', ...
        systemLabelLocal(sysName)));

    legend(ax, ...
        {'Reference min', 'Reference max', ...
         'Coordinated min', 'Coordinated max'}, ...
        'Location', 'best');

    ylim(ax, [ ...
        min([refMin; coordMin; cfg.Vmin]) - 0.02, ...
        max([refMax; coordMax; cfg.Vmax]) + 0.02]);

    removeAxesToolbarLocal(ax);

    cleanStem = sprintf('fig_%s_daily_voltage_profiles', sysName);

    % Former filename retained only for backward compatibility with the
    % current LaTeX source. It is not shown inside the paper.
    legacyStem = sprintf('fig_%s_safety_gated_voltage_envelope', sysName);

    saveFigureCopiesLocal(fig, cfg.outputDir, cleanStem, legacyStem);
    close(fig);
end

function makeViolationFigureLocal(sysName, Vref, Vcoord, cfg)
    t = 1:size(Vref, 1);

    refViol = sum(Vref < cfg.Vmin | Vref > cfg.Vmax, 2);
    coordViol = sum(Vcoord < cfg.Vmin | Vcoord > cfg.Vmax, 2);

    fig = figure( ...
        'Color', 'w', ...
        'Position', [100 100 900 600], ...
        'MenuBar', 'none', ...
        'ToolBar', 'none');

    ax = axes(fig);
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');

    plot(ax, t, refViol, '--', 'LineWidth', 1.8);
    plot(ax, t, coordViol, '-', 'LineWidth', 2.2);

    xlabel(ax, 'Time step');
    ylabel(ax, 'Number of violating buses');

    title(ax, sprintf('%s voltage-limit violations', ...
        systemLabelLocal(sysName)));

    legend(ax, ...
        {'Reference operation', 'Coordinated strategy'}, ...
        'Location', 'best');

    ylim(ax, [0, max([refViol; coordViol]) + 2]);

    removeAxesToolbarLocal(ax);

    cleanStem = sprintf('fig_%s_voltage_limit_violations', sysName);

    % Former filename retained only for backward compatibility with the
    % current LaTeX source. It is not shown inside the paper.
    legacyStem = ...
        sprintf('fig_%s_safety_gated_violation_count_per_step', sysName);

    saveFigureCopiesLocal(fig, cfg.outputDir, cleanStem, legacyStem);
    close(fig);
end

function removeAxesToolbarLocal(ax)
    % Prevent interactive MATLAB toolbar icons from appearing in exports.
    try
        disableDefaultInteractivity(ax);
    catch
        % Continue safely if this function is unavailable.
    end

    try
        axtoolbar(ax, {});
    catch
        % Continue safely if this function is unavailable.
    end
end

function saveFigureCopiesLocal(fig, outputDir, cleanStem, legacyStem)
    cleanBase = fullfile(outputDir, cleanStem);
    legacyBase = fullfile(outputDir, legacyStem);

    drawnow;

    savefig(fig, [cleanBase '.fig']);
    exportgraphics(fig, [cleanBase '.png'], 'Resolution', 300);

    % Copy the clean files to the former names so existing LaTeX references
    % continue to work until the manuscript file is updated.
    copyfile([cleanBase '.fig'], [legacyBase '.fig']);
    copyfile([cleanBase '.png'], [legacyBase '.png']);

    fprintf('Saved:\n  %s.png\n', cleanBase);
end

function label = systemLabelLocal(sysName)
    switch lower(sysName)
        case '33bus'
            label = '33-bus';
        case '69bus'
            label = '69-bus';
        otherwise
            label = sysName;
    end
end

function writeLatexSummaryLocal(T, outFile)
    fid = fopen(outFile, 'w');

    if fid < 0
        error('Cannot write %s', outFile);
    end

    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '\\begin{table}[!t]\n');
    fprintf(fid, ['\\caption{Daily voltage-quality comparison for the ', ...
                  'modified 33-bus and 69-bus systems.}\n']);
    fprintf(fid, '\\label{tab:daily_voltage_quality_summary}\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{llccccc}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, ...
        'System & Case & Min. V & Max. V & Mean V & $J_V$ & Viol. \\\\ \n');
    fprintf(fid, '\\hline\n');

    for i = 1:height(T)
        fprintf(fid, '%s & %s & %.4f & %.4f & %.4f & %.4f & %d \\\\ \n', ...
            T.System{i}, T.Case{i}, T.MinV(i), T.MaxV(i), T.MeanV(i), ...
            T.TotalVoltageDeviation(i), T.TotalViolationCount(i));
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
end
