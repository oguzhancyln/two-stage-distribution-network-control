%% find_inverter_apparent_power_setting.m
% Searches MATLAB source files for the inverter apparent-power rating
% definition used in the hybrid-voltage-control implementation.
%
% Run this script from the main project folder. Copy the printed matches
% and inspect the line that defines the inverter rating (S_k / Sinv / etc.).

clear; clc;

projectRoot = pwd;
files = dir(fullfile(projectRoot, '**', '*.m'));

patterns = [ ...
    "Sinv", ...
    "S_inv", ...
    "Smax", ...
    "S_max", ...
    "Srated", ...
    "S_rated", ...
    "apparent", ...
    "rating", ...
    "qmax", ...
    "Qmax", ...
    "pfmin", ...
    "pf_min", ...
    "pvRated", ...
    "PVrated", ...
    "P_rated", ...
    "Ppv" ...
];

fprintf('Searching MATLAB files under:\n  %s\n\n', projectRoot);

matchCount = 0;

for f = 1:numel(files)
    filePath = fullfile(files(f).folder, files(f).name);

    try
        txt = fileread(filePath);
    catch
        continue;
    end

    lines = splitlines(string(txt));
    printedHeader = false;

    for i = 1:numel(lines)
        lineText = lines(i);

        if any(contains(lineText, patterns, 'IgnoreCase', true))
            if ~printedHeader
                fprintf('\n============================================================\n');
                fprintf('%s\n', filePath);
                fprintf('============================================================\n');
                printedHeader = true;
            end

            fprintf('%5d | %s\n', i, lineText);
            matchCount = matchCount + 1;
        end
    end
end

fprintf('\n============================================================\n');
fprintf('Total matching lines: %d\n', matchCount);
fprintf(['Inspect the lines that define the inverter apparent-power ', ...
         'rating and the dynamic Qmax expression.\n']);
