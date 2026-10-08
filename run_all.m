% RUN_ALL  Generate manuscript Figures 2-8 and A1.
%
% Runs the eight supplied figure scripts in manuscript order. Each script
% uses its own parameters and resolves data paths relative to its location.
% Figures 3 and 6 include time-dependent numerical calculations; Figure 3
% also records an animation using a graphical MATLAB session.
%
% Run from the extracted project folder with the command run_all. PNG files
% and the GIF are written to figures/. The session log records the local
% MATLAB version, installed products and command-window output. Each figure
% script retains its original workspace-clearing behavior.
%
% Input: matlab/FIG_2.m through FIG_8.m, matlab/FIG_A1.m and data/.
% Output: eight PNG files, DEVOL_2Fronts.gif and run_log.txt in figures/.
% Figure 1 has no supplied source.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

if ~exist(fullfile(fileparts(mfilename('fullpath')), 'figures'), 'dir')
    mkdir(fullfile(fileparts(mfilename('fullpath')), 'figures'));
end
diary off
diary(fullfile(fileparts(mfilename('fullpath')), 'figures', 'run_log.txt'));
fprintf('Run started: %s\n', datestr(now));
ver

% Separate calls are needed because the figure scripts clear their workspace.
fprintf('\nFigure 2\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_2.m'));
fprintf('\nFigure 3\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_3.m'));
fprintf('\nFigure 4\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_4.m'));
fprintf('\nFigure 5\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_5.m'));
fprintf('\nFigure 6\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_6.m'));
fprintf('\nFigure 7\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_7.m'));
fprintf('\nFigure 8\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_8.m'));
fprintf('\nFigure A1\n');
run(fullfile(fileparts(mfilename('fullpath')), 'matlab', 'FIG_A1.m'));

fprintf('\nRun finished: %s\n', datestr(now));
fprintf('Figures saved in %s\n', fullfile(fileparts(mfilename('fullpath')), 'figures'));
diary off
