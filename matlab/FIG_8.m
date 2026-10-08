% FIG_8  Permeability estimates for the gypsum dehydration experiment.
%
% Applies the manuscript's single- and two-front solutions to the gypsum
% dehydration experiment of Fusseis et al. (2012). Panel (a) compares
% permeabilities inferred from hydraulic diffusion, a single dehydration
% or vaporization front, and two coupled dehydration/vaporization fronts.
%
% The one-front estimates use the measured coefficient 8.29e-11 m^2/s.
% The two-front estimates search a logarithmic permeability grid for a
% front position of 1 mm at 12000 s. Panel (b) plots digitized porosity
% profiles at 600, 4200 and 7800 s, redrawn after Fusseis et al. (2012).
%
% Input: data/DATA_Fusseis_plots.mat; x1,y1, x2,y2 and x3,y3 contain
%        distance in micrometres and porosity in percent.
% Output: figures/FIG_8.png.
%
% Reference:
% Fusseis, F., Schrank, C., Liu, J., Karrech, A., Llana-Funez, S.,
% Xiao, X., and Regenauer-Lieb, K. (2012). Pore formation during dehydration
% of a polycrystalline gypsum sample observed and quantified in a time-series
% synchrotron X-ray micro-tomography experiment. Solid Earth, 3, 71-86.
% https://doi.org/10.5194/se-3-71-2012
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear; close all; clc;

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FS_AXIS   = 16;   % axes labels and ticks
FS_TITLE  = 16;   % titles
FS_TEXT   = 16;   % annotations on main plot
FS_LEG    = 16;   % custom legend table text
FS_INSET  = 16;   % inset fonts

LW_MAIN   = 3.0;
LW_INSET  = 1.8;
LW_BOX    = 1.0;
LW_DIV    = 0.9;

DP_vec          = (10:100)*1e6;    % pressure difference, Pa
Diffusivity     = 8.29e-11;        % measured diffusivity in Fusseis experiment
Drhot           = 100;             % density jump for dehydration front
Drhotv          = 320;             % density jump for vapor front
DPv             = 0.14e6;          % pressure difference for water-vapor transition, Pa
rhof            = 1e3;             % water density, kg/m^3
rhof_v          = 0.52;            % vapor density, kg/m^3
viscosity       = 3e-4;            % water viscosity, Pa s
viscosity_v     = 1.25e-5;         % vapor viscosity, Pa s
kin_viscosity_v = 2.4e-5;          % vapor kinematic viscosity, m^2/s
K               = 50e9;            % poroelastic modulus, Pa
t_exp           = 12e3;            % experimental time, s

k_hyd = zeros(size(DP_vec));
k_eff = zeros(size(DP_vec));
k_v1  = zeros(size(DP_vec));
kf1   = zeros(size(DP_vec));
kf2   = zeros(size(DP_vec));

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', ...
    'DATA_Fusseis_plots.mat'), 'x1', 'y1', 'x2', 'y2', 'x3', 'y3');
x_exp1 = x1;  y_exp1 = y1;
x_exp2 = x2;  y_exp2 = y2;
x_exp3 = x3;  y_exp3 = y3;
clear x1 y1 x2 y2 x3 y3

% one-front estimates
for it = 1:length(DP_vec)
    DP = DP_vec(it);
    k_hyd(it) = Diffusivity/K*viscosity;
    k_eff(it) = Diffusivity/(2*DP)/rhof * viscosity * Drhot;
    k_v1(it)  = Diffusivity/(2*DPv) * kin_viscosity_v * Drhotv;
end

% two-front estimates
for it = 1:length(DP_vec)
    DP  = DP_vec(it);
    k2f = 10.^linspace(-25,-13,1e4);

    Deff2 = rhof_v/Drhotv .* k2f/viscosity_v * DPv;
    Deff1 = rhof  /Drhot  .* k2f/viscosity   * DP;

    gam     = Drhot/Drhotv;
    r       = (-(Deff2 + gam*Deff1) + sqrt((Deff2 + gam*Deff1).^2 + 4*Deff1.*Deff2)) ./ (2*Deff1);

    Lam1 = sqrt(2*Deff1./(1-r));
    Lam2 = r.*Lam1;

    xfront1 = Lam1.*sqrt(t_exp);
    xfront2 = Lam2.*sqrt(t_exp);

    [~,if1] = min(abs(xfront1 - 1e-3));
    kf1(it) = k2f(if1);

    [~,if2] = min(abs(xfront2 - 1e-3));
    kf2(it) = k2f(if2);
end

fig = figure(1); clf;
set(fig,'Color','w','Position',[100 100 1280 700]);

ax_main = axes('Position',[0.08 0.12 0.5 0.79]);
hold(ax_main,'on'); box(ax_main,'on');

p1 = plot(ax_main, DP_vec/1e6, k_hyd,      '-k',  'LineWidth',LW_MAIN);
p2 = plot(ax_main, DP_vec/1e6, 0.9*k_eff, '--b', 'LineWidth',LW_MAIN);
p3 = plot(ax_main, DP_vec/1e6, k_v1,       '--r', 'LineWidth',LW_MAIN);
p4 = plot(ax_main, DP_vec/1e6, kf2,        '-r',  'LineWidth',LW_MAIN);
p5 = plot(ax_main, DP_vec/1e6, kf1,        '-b',  'LineWidth',LW_MAIN);

set(ax_main,'XScale','log','YScale','log');
set(ax_main,'YTick',[1e-25 1e-23 1e-21 1e-19 1e-17 1e-15 1e-13 1e-11]);
axis(ax_main,[min(DP_vec)/1e6 max(DP_vec)/1e6 1e-25 1e-13]);

grid(ax_main,'on');
set(ax_main,'FontSize',FS_AXIS*1.1,'LineWidth',1.0,'Layer','top');

xlabel(ax_main,'Fluid pressure difference (MPa)','FontSize',FS_AXIS*1.1);
ylabel(ax_main,'Permeability (m$^2$)','FontSize',FS_AXIS*1.1);
title(ax_main, {'$\textbf{(a) Permeability estimates for different analytical solutions}$'},'Interpreter','latex','FontSize',FS_TITLE*1.1);

text(ax_main, 11, 1e-24, '\textbf{Non-reactive hydraulic diffusion front}', ...
    'Interpreter','latex', 'FontSize',FS_TEXT, 'Color','k');

text(ax_main, 11, 2.5e-23, '\textbf{Single dehydration front}', ...
    'Interpreter','latex', 'FontSize',FS_TEXT, 'Color','b');

text(ax_main, 11, 1e-18, '\textbf{Single vaporization front}', ...
    'Interpreter','latex', 'FontSize',FS_TEXT, 'Color','r');

text(ax_main, 32, 8e-23, '\textbf{Two fronts: leading dehydration front}', ...
    'Interpreter','latex', 'FontSize',FS_TEXT, 'Color','b');

text(ax_main, 32, 6e-15, '\textbf{Two fronts: trailing vaporization front}', ...
    'Interpreter','latex', 'FontSize',FS_TEXT, 'Color','r');

% Legend
ax_leg = axes('Position',[0.67-0.08 0.06 0.32 0.32]);
axis(ax_leg,[0 1 0 1]);
axis(ax_leg,'off');
hold(ax_leg,'on');

rectangle(ax_leg,'Position',[0 0.2 1 0.8], ...
    'EdgeColor','k', 'LineWidth',LW_BOX, 'FaceColor','w');

line(ax_leg,[0.22 0.22],[0.2 1], 'Color','k', 'LineWidth',LW_DIV);
line(ax_leg,[0.55 0.55],[0.2 1], 'Color','k', 'LineWidth',LW_DIV);
line(ax_leg,[0 1],[0.80 0.80], 'Color','k', 'LineWidth',LW_DIV);
line(ax_leg,[0 1],[0.67 0.67], 'Color','k', 'LineWidth',LW_DIV);
line(ax_leg,[0 1],[0.425 0.425], 'Color','k', 'LineWidth',LW_DIV);

text(ax_leg,0.38,1.05,'\textbf{Legend and assumptions of (a)}', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG*1.1);

text(ax_leg,0.385,0.92,'\textbf{Measured}', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.385,0.85,'\textbf{front linked to}', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

text(ax_leg,0.775,0.92,'\textbf{Analytical}', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,0.85,'\textbf{solution}', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

y = [0.725 0.60 0.50 0.35 0.25];

plot(ax_leg,[0.05 0.17],[y(1) y(1)],'-k','LineWidth',LW_INSET);
text(ax_leg,0.385,y(1),'dehydration', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,y(1),'hydraulic diffusion', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

plot(ax_leg,[0.05 0.17],[y(2) y(2)],'--','Color','b','LineWidth',LW_INSET);
text(ax_leg,0.385,y(2),'dehydration', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,y(2),'single', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

plot(ax_leg,[0.05 0.17],[y(3) y(3)],'--','Color','r','LineWidth',LW_INSET);
text(ax_leg,0.385,y(3),'vaporization', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,y(3),'reaction front', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

plot(ax_leg,[0.05 0.17],[y(4) y(4)],'-','Color','b','LineWidth',LW_INSET);
text(ax_leg,0.385,y(4),'dehydration', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,y(4),'two coupled', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

plot(ax_leg,[0.05 0.17],[y(5) y(5)],'-','Color','r','LineWidth',LW_INSET);
text(ax_leg,0.385,y(5),'vaporization', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);
text(ax_leg,0.775,y(5),'reaction fronts', ...
    'Interpreter','latex', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontSize',FS_LEG);

% Measured porosity profiles
ax_inset = axes('Position',[0.71-0.08 0.54 0.28 0.37]);
hold(ax_inset,'on'); box(ax_inset,'on');

plot(ax_inset, x_exp1, y_exp1,   'k--', 'LineWidth', LW_INSET);
plot(ax_inset, x_exp2, y_exp2,   'k-.', 'LineWidth', LW_INSET);
plot(ax_inset, x_exp3, y_exp3-1, 'k-',  'LineWidth', LW_INSET);

xlabel(ax_inset,'Distance ($\mu$m)','Interpreter','latex','FontSize',FS_INSET);
ylabel(ax_inset,'Porosity (\%)','Interpreter','latex','FontSize',FS_INSET);
title(ax_inset,'\textbf{(b) Measured porosity profiles}','Interpreter','latex','FontSize',FS_INSET);

axis(ax_inset,[0 570 0 34]);
grid(ax_inset,'on');
set(ax_inset,'FontSize',FS_INSET,'LineWidth',1.0,'Layer','top');

legend(ax_inset,{'600 s','4200 s','7800 s'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',FS_INSET-2);

drawnow;
saveas(fig,fullfile(fig_dir,'FIG_8.png'));
