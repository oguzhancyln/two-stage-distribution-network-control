%% regenerate_figures_09_10_network_control_convergence.m
% Regenerates Figures 9 and 10 from the saved optimizer-benchmark CSV.
% No optimization is rerun.
%
% Required input:
%   test03_benchmark_convergence.csv
%
% Expected CSV columns:
%   System, Method, Run, Window, Iteration, BestObjective
%
% The script uses Window = 3, corresponding to time steps 73--96.
% Output PNG files are written under the figures folder with the same
% filenames already referenced by the LaTeX manuscript.

clear; clc; close all;

%% USER SETTING
% Run this script from the main paper/project folder.
projectRoot = pwd;

% Output folder used by the LaTeX manuscript.
figuresDir = fullfile(projectRoot, 'figures');
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

%% FIND THE SAVED CONVERGENCE CSV
csvName = 'test03_benchmark_convergence.csv';

% First check the current project folder.
csvFile = fullfile(projectRoot, csvName);

% If the CSV is not in the project root, search subfolders recursively.
if ~isfile(csvFile)
    matches = dir(fullfile(projectRoot, '**', csvName));
    if isempty(matches)
        error(['Could not find ', csvName, ...
            '. Copy the CSV into the project folder or a subfolder, ', ...
            'then run the script again.']);
    end

    % Use the most recently modified matching file if more than one exists.
    [~, idxNewest] = max([matches.datenum]);
    csvFile = fullfile(matches(idxNewest).folder, matches(idxNewest).name);
end

fprintf('Reading convergence data from:\n  %s\n\n', csvFile);

%% LOAD AND VALIDATE DATA
T = readtable(csvFile, 'TextType', 'string');

requiredVars = ["System", "Method", "Run", "Window", "Iteration", "BestObjective"];
missingVars = requiredVars(~ismember(requiredVars, string(T.Properties.VariableNames)));

if ~isempty(missingVars)
    error('Missing required CSV column(s): %s', strjoin(missingVars, ', '));
end

T.System = strtrim(string(T.System));
T.Method = upper(strtrim(string(T.Method)));

% Window 3 is the late-day network-control interval: time steps 73--96.
lateDayWindow = 3;

systems = ["33bus", "69bus"];
methods = ["GWO", "PSO", "DE"];

outputFiles = [ ...
    "fig_test03_33bus_late_day_dnrcap_convergence.png", ...
    "fig_test03_69bus_late_day_dnrcap_convergence_correct.png" ...
];

figureTitles = [ ...
    "33-bus: late-day network-control convergence", ...
    "69-bus: late-day network-control convergence" ...
];

%% GENERATE FIGURES 9 AND 10
for s = 1:numel(systems)
    systemName = systems(s);

    fig = figure( ...
        'Color', 'w', ...
        'Position', [100, 100, 980, 610]);

    ax = axes(fig);
    hold(ax, 'on');

    for m = 1:numel(methods)
        methodName = methods(m);

        rows = T.System == systemName & ...
               T.Method == methodName & ...
               T.Window == lateDayWindow;

        Tm = T(rows, :);

        if isempty(Tm)
            error('No data found for %s / %s / Window %d.', ...
                systemName, methodName, lateDayWindow);
        end

        iterations = unique(Tm.Iteration);
        meanBest = zeros(size(iterations));

        for k = 1:numel(iterations)
            iterRows = Tm.Iteration == iterations(k);
            meanBest(k) = mean(Tm.BestObjective(iterRows), 'omitnan');
        end

        semilogy(ax, iterations, meanBest, 'LineWidth', 2.3);

        fprintf('%s / %s: %d runs, %d iterations\n', ...
            systemName, methodName, ...
            numel(unique(Tm.Run)), numel(iterations));
    end

    grid(ax, 'on');
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'off';
    ax.FontSize = 12;
    ax.LineWidth = 1.0;

    xlabel(ax, 'Iteration', 'FontSize', 13);
    ylabel(ax, 'Mean best-so-far objective', 'FontSize', 13);
    title(ax, figureTitles(s), 'FontSize', 15, 'FontWeight', 'bold');
    legend(ax, methods, 'Location', 'northeast', 'FontSize', 12);

    xlim(ax, [1, max(T.Iteration)]);

    pngFile = fullfile(figuresDir, outputFiles(s));
    exportgraphics(fig, pngFile, 'Resolution', 300);

    fprintf('Saved:\n  %s\n\n', pngFile);

    close(fig);
end

fprintf('Done. Figures 9 and 10 were regenerated without rerunning optimization.\n');
