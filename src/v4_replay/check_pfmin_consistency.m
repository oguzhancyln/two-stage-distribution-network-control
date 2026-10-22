%% check_pfmin_consistency.m
% Scans MATLAB files under the current project folder and reports every
% cfg.PFmin / cfg.pfmin assignment. It does NOT modify any file.

clear; clc;

projectRoot = pwd;
files = dir(fullfile(projectRoot, '**', '*.m'));

priorityFiles = lower({ ...
    'runtime_benchmark.m', ...
    'runtime_benchmark_v2.m', ...
    'stage4_losses_switching.m' ...
});

fprintf('Scanning MATLAB files under:\n  %s\n\n', projectRoot);

matchCount = 0;
priorityCount = 0;

for k = 1:numel(files)
    filePath = fullfile(files(k).folder, files(k).name);

    try
        txt = fileread(filePath);
    catch
        continue;
    end

    lines = splitlines(string(txt));
    printedHeader = false;

    for i = 1:numel(lines)
        lineText = lines(i);

        % Search common assignment forms:
        % cfg.PFmin = 0.90;
        % cfg.pfmin = 0.95;
        % PFmin = 0.90;
        % pfmin = 0.95;
        if ~isempty(regexp(char(lineText), ...
                '(?i)\b(?:cfg\.)?pfmin\s*=\s*[0-9.]+', 'once'))

            if ~printedHeader
                fprintf('\n============================================================\n');
                fprintf('%s\n', filePath);
                fprintf('============================================================\n');
                printedHeader = true;
            end

            isPriority = any(strcmpi(files(k).name, priorityFiles));

            if isPriority
                fprintf('PRIORITY ');
                priorityCount = priorityCount + 1;
            else
                fprintf('         ');
            end

            fprintf('%5d | %s\n', i, strtrim(lineText));
            matchCount = matchCount + 1;
        end
    end
end

fprintf('\n============================================================\n');
fprintf('Total PFmin assignments found: %d\n', matchCount);
fprintf('Priority-file matches found:   %d\n', priorityCount);
fprintf(['\nFor the final paper helper analyses, confirm that the active ', ...
         'assignments in runtime_benchmark*.m and stage4_losses_switching.m ', ...
         'use 0.90.\n']);
