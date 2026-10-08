% FIG_A1  Two-front geometry and the pressure-density relation.
%
% Draws the schematic geometry used in Appendix A. The pressure profile
% contains a leading front x1 and a trailing front x2, with
% P_LB < P_F2 < P_F1 = P_ini.
%
% The second panel shows the piecewise-linear constitutive relation with
% two density jumps at fixed pressure. Generic variables u and A represent
% pressure and total density in the labels. Front locations and material
% parameters are illustrative; this script does not solve a time-dependent
% transport problem.
%
% Input: schematic parameters specified in this script.
% Output: figures/FIG_A1.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear; clc; close all;

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

c1 = 5;
c2 = 8;
c3 = 10;

uF1 = 1.2;
uF2 = 2.6;

AF1 = 3.0;

DeltaA1 = 3.0;
DeltaA2 = 4.0;

% Continuity along the middle branch
AF2 = (AF1 + DeltaA1) + c2*(uF2 - uF1);

% Single-phase branches
u1 = linspace(0.0, uF1, 200);
u2 = linspace(uF1, uF2, 200);
u3 = linspace(uF2, 3.4, 200);

% Piecewise-linear A(u)
A1 = AF1 + c1*(u1 - uF1);
A2 = (AF1 + DeltaA1) + c2*(u2 - uF1);
A3 = (AF2 + DeltaA2) + c3*(u3 - uF2);

ulb  = 0.8;
uini = uF2*1.0000000001;

figure('Color','w');

yticks_vals   = [ulb, uF1, uF2];
yticks_labels = {'$P_{LB}$','$P_{F2}$','$P_{F1}=P_{ini}$'};

subplot(121)
hold on; box on;

plot([0,0.3],[ulb, uF1], 'b-', 'LineWidth', 2.5)
plot([0.3,0.5],[uF1, uF2], 'r-', 'LineWidth', 2.5)
plot([0.5,1],[uini, uini], 'r--', 'LineWidth', 2.5)

xline(0.3,  'b--', 'LineWidth', 1.2);
xline(0.5,  'r--', 'LineWidth', 1.2);
text(0.32, uF2 + 0.4, '$x_2(t)$', ...
    'Interpreter','latex', 'FontSize', 20,'Color','b');
text(0.52, uF2 + 0.4, '$x_1(t)$', ...
    'Interpreter','latex', 'FontSize', 20,'Color','r');

yline(uini, 'g:', 'LineWidth', 1.2);
yline(ulb,  'b:', 'LineWidth', 1.2);
yline(uF1,  'k:', 'LineWidth', 1.2);
yline(uF2,  'k:', 'LineWidth', 1.2);

legend({'Region 2', 'Region 1', 'Region 0', ...
        'Front 2', 'Front 1'}, ...
        'Interpreter','latex', 'Location','southeast')

ylim([0, max(u3)*1.05]);
xlim([0, 1]);

grid on;
set(gca, 'FontSize', 18, 'LineWidth', 1.0);

set(gca, ...
    'YTick', yticks_vals, ...
    'YTickLabel', yticks_labels, ...
    'TickLabelInterpreter', 'latex', ...
    'XTickLabel', []);

ylabel('$P$', 'Interpreter','latex', 'FontSize', 18);
xlabel('$x$', 'Interpreter','latex', 'FontSize', 18);

title('(a) Fluid pressure versus spatial position', ...
    'Interpreter','latex', 'FontSize', 18);

subplot(122)
hold on; box on;

% Single-phase branches
plot(A1, u1, 'b-', 'LineWidth', 2.5);
plot(A2, u2, 'r-', 'LineWidth', 2.5);
plot(A3, u3, 'g-', 'LineWidth', 2.5);

% Density jumps at fixed pressure
plot([AF1, AF1 + DeltaA1], [uF1, uF1], 'k:', 'LineWidth', 2.5);
plot([AF2, AF2 + DeltaA2], [uF2, uF2], ' :', 'LineWidth', 2.5,'Color',[1 1 1]*0.5);

Alb  = interp1(u1, A1, ulb);
Aini = interp1(u3, A3, uini);

plot(Aini, uini, 'go', 'MarkerFaceColor','g', 'MarkerSize', 7);
plot(Alb,  ulb,  'bo', 'MarkerFaceColor','b', 'MarkerSize', 7);

yline(uF1,  'k:', 'LineWidth', 1.2);
yline(uF2,  'r:', 'LineWidth', 1.2);
yline(uini, 'r:', 'LineWidth', 1.2);
yline(ulb,  'b:', 'LineWidth', 1.2);

xlabel('$\rho_{total}(P)$', 'Interpreter','latex', 'FontSize', 18);
ylabel('$P$', 'Interpreter','latex', 'FontSize', 18);

legend({'Assemblage 2', 'Assemblage 1', 'Assemblage 0', ...
        'Reaction 2', 'Reaction 1'}, ...
        'Interpreter','latex', 'Location','southeast');

text(AF1 + 0.0*DeltaA1, uF1 - 0.08, '$\Delta \rho_{total,2}$', ...
    'Interpreter','latex', 'FontSize', 20);
text(AF2 + 0.0*DeltaA2, uF2 - 0.08, '$\Delta \rho_{total,1}$', ...
    'Interpreter','latex', 'FontSize', 20);

grid on;
set(gca, 'FontSize', 18, 'LineWidth', 1.0);

xlim([0, max(A3)*1.05]);
ylim([0, max(u3)*1.05]);

set(gca, ...
    'YTick', yticks_vals, ...
    'YTickLabel', yticks_labels, ...
    'TickLabelInterpreter', 'latex', ...
    'XTickLabel', []);

title('(b) Fluid pressure versus total density', ...
    'Interpreter','latex', 'FontSize', 18);

set(gcf,'Position', [169 77 1248 677],'color',[1 1 1])
saveas(gcf, fullfile(fig_dir, 'FIG_A1.png'));
