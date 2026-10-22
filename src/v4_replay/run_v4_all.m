%% run_v4_all.m
% Convenience driver with MATLAB-safe short filenames.
clear; clc;

thisFile = mfilename('fullpath');
if isempty(thisFile)
    scriptDir = pwd;
else
    scriptDir = fileparts(thisFile);
end

addpath(scriptDir);
cd(scriptDir);

fprintf('\n[1/2] Running V4 consistent-Qmax ablation...\n');
run(fullfile(scriptDir,'v4_ablation.m'));

fprintf('\n[2/2] Generating V4 figures and summary...\n');
run(fullfile(scriptDir,'v4_figures.m'));

fprintf('\nDone. Outputs are under:\n');
fprintf('  %s\n', fullfile(scriptDir,'results_ablation_33_69_v4_consistent_qmax'));
fprintf('  %s\n', fullfile(scriptDir,'results_test01c_consistent_qmax_envelopes'));
