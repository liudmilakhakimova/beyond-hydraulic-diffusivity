% FIG_6  Hydraulic diffusion with zero, one and two reactions.
%
% Solves the same one-dimensional mass-conservation and Darcy-flow problem
% as FIG_3 for three pressure-density relations: a continuous linear law,
% a law with one reaction interval, and a law with two reaction intervals.
%
% The left column shows pressure profiles and propagation-distance markers;
% the right column shows the corresponding constitutive relations. The
% comparison illustrates reaction-induced retardation of the pressure
% perturbation.
%
% Input: physical and numerical parameters specified in this script.
% Output: figures/FIG_6.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear; close all; clc;

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure(1); clf;
set(figure(1),'Color','w','Position',[10 10 1380 1200]);

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FS = 16;

LW_PROFILE = 2.8;
LW_FRONT   = 2.8;
LW_LUT     = 2.8;

figure(1); clf; colormap jet

for icase = 1:3

    % Physical parameters
    k_etaf = 1;        % permeability / fluid viscosity, m^2 / Pa / s
    Lx     = 1;        % model length, m
    rhof   = 1;        % fluid density, kg / m^3
    t_pert = 1e-1;     % experiment time, s

    % Constitutive relation
    if icase == 1
        rhotp = [4.0 4.0];   % total density after the 1st and 2nd reactions
        rhotm = [4.0 4.0];   % total density before the 1st and 2nd reactions
    elseif icase == 2
        rhotp = [2.5 4.0];
        rhotm = [2.5 2.5];
    elseif icase == 3
        rhotp = [2.0 4.0];
        rhotm = [1.0 2.5];
    end

    beta      = [1 1 1];   % compressibility in the 1st, 2nd, and 3rd assemblages
    Pf_react2 = 1.5;
    Pf_react1 = Pf_react2 - (rhotm(2) - rhotp(1)) / (rhotm(2) * beta(2));

    % Piecewise-linear pressure-density relation
    rhot_min = 0.1 * rhotm(1);
    rhot_max = 1.1 * rhotp(2);
    lutN     = 40000;

    rhot_lut = linspace(rhot_min, rhot_max, lutN);
    Pf_lut   = zeros(size(rhot_lut));

    dPdrho = [1/(rhotm(1)*beta(1)), 1/(rhotm(2)*beta(2)), 1/(rhotp(2)*beta(3))];

    for i = 1:length(rhot_lut)
        rho = rhot_lut(i);

        if rho <= rhotm(1)
            % first assemblage
            Pf_lut(i) = Pf_react1 + (rho - rhotm(1)) * dPdrho(1);
        elseif rho <= rhotp(1)
            % first reaction
            Pf_lut(i) = Pf_react1;
        elseif rho <= rhotm(2)
            % second assemblage
            Pf_lut(i) = Pf_react1 + (rho - rhotp(1)) * dPdrho(2);
        elseif rho <= rhotp(2)
            % second reaction
            Pf_lut(i) = Pf_react2;
        else
            % third assemblage
            Pf_lut(i) = Pf_react2 + (rho - rhotp(2)) * dPdrho(3);
        end
    end

    % Numerical parameters
    nx   = 250;
    nt   = 65000;
    nout = 1000;

    dx = Lx/(nx-1);
    x  = linspace(0, Lx, nx);
    dt = dx^2 / (rhof * k_etaf * max(dPdrho)) / 2.1;

    % Initial and boundary conditions
    rhot_lb  = 0.5;               % total density at the left boundary
    rhot_ini = 1.0 * rhotp(2);      % total density of the initial assemblage
    Pf_lb    = interp1(rhot_lut, Pf_lut, rhot_lb, 'linear', 'extrap');

    rhot    = rhot_ini * ones(size(x));
    rhot(1) = rhot_lb;

    if icase == 3
        Pf_lb    = 0.7;
        rhot_lb  = interp1(Pf_lut(1:100), rhot_lut(1:100), Pf_lb, 'linear', 'extrap');
        rhot(1)  = rhot_lb;
    end

    % Quasi-stationary two-front solution
    if icase == 1
        DeltaRhot = [rhotp(1)-rhotm(1), rhotp(2)-rhotm(2)];
        DeltaPf   = [Pf_react1-Pf_lb, Pf_react2-Pf_react1];

        Liu1  = DeltaRhot(1)/(rhotm(1)*beta(1)*DeltaPf(1));
        Liu2  = DeltaRhot(2)/(rhotm(2)*beta(2)*DeltaPf(2));
        Deff1 = rhof*k_etaf/(rhotm(1)*beta(1)*Liu1);
        Deff2 = rhof*k_etaf/(rhotm(2)*beta(2)*Liu2);
    else
        DeltaRhot = [rhotp(2)-rhotm(2), rhotp(1)-rhotm(1)];
        DeltaPf   = [abs(Pf_react2-Pf_lb), abs(Pf_react2-Pf_react1)];

        Liu1  = DeltaRhot(1)/(rhotp(2)*beta(3)*DeltaPf(1));
        Liu2  = DeltaRhot(2)/(rhotp(1)*beta(2)*DeltaPf(2));
        Deff1 = rhof*k_etaf/(rhotp(2)*beta(3)*Liu1);
        Deff2 = rhof*k_etaf/(rhotp(1)*beta(2)*Liu2);
    end

    beta_qs = DeltaRhot(2)/DeltaRhot(1);
    r       = (-(Deff1 + beta_qs*Deff2) + sqrt((Deff1 + beta_qs*Deff2)^2 + 4*Deff1*Deff2)) / (2*Deff2);
    lambda2 = sqrt(2*Deff2/(1-r));
    lambda1 = r*lambda2;

    time     = [];
    xf01_num = [];
    xf01_ana = [];
    xf1_num  = [];
    xf2_num  = [];
    xf1_ana  = [];
    xf2_ana  = [];

    it  = 0;
    Pf0 = interp1(rhot_lut, Pf_lut, rhot, 'linear', 'extrap');

    while it*dt < t_pert
        it = it + 1;

        % pressure from LUT
        Pf = interp1(rhot_lut, Pf_lut, rhot, 'linear', 'extrap');

        % Darcy flux + total mass conservation
        qD = -k_etaf * diff(Pf) / dx;
        rhot(2:end-1) = rhot(2:end-1) - dt * diff(rhof*qD) / dx;

        % Reaction-front positions
        if mod(it, nout) == 0
            time(end+1) = it*dt;

            % numerical
            if icase == 2
                idx1 = find(Pf <= Pf_react1, 1, 'last');
                xf01_num(end+1) = interp1(Pf(idx1:idx1+1), x(idx1:idx1+1), Pf_react1, 'linear', 'extrap');
                xf01_num_last   = xf01_num(end);
            elseif icase == 3
                idx1 = find(Pf <= Pf_react1, 1, 'last');
                xf1_num(end+1) = interp1(Pf(idx1:idx1+1), x(idx1:idx1+1), Pf_react1, 'linear', 'extrap');

                idx2 = find(Pf <= Pf_react2, 1, 'last');
                xf2_num(end+1) = interp1(Pf(idx2:idx2+1), x(idx2:idx2+1), Pf_react2, 'linear', 'extrap');
            end

            % analytical
            xf1_ana(end+1)  = sqrt(2*k_etaf*time(end));
            if icase==1; xf1_ana_last    = xf1_ana(end); end
            if icase==2
                xf01_ana(end+1) = sqrt(2*k_etaf*rhof*(Pf_react1-Pf_lb)/DeltaRhot(1)) * sqrt(time(end));
                dP_case2        = Pf_react1-Pf_lb;
                Drho_case2      = DeltaRhot(1);
            elseif icase==3
                xf01_ana(end+1) = sqrt(2*k_etaf*rhof*(dP_case2)/Drho_case2) * sqrt(time(end));
            end

            col_f0_num  = [0.20 0.70 0.80];   % cyan
            col_f01_num = [0.50 0.50 0.50];   % grey
            col_f1_num  = [0.00 0.35 0.70];   % dark blue
            col_f2_num  = [0.75 0.15 0.10];   % dark red

            if icase == 1
                figure(1); clf;
            else
                figure(1);
            end
            set(gcf,'Color','w');

            % Constitutive relation
            subplot(3,2,2 + (icase-1)*2); hold on; box on;
            plot(rhot_lut, -(Pf_lut-Pf_react2)/(Pf_lb-Pf_react2), 'k-', 'LineWidth', LW_LUT, 'HandleVisibility', 'off');
            plot(rhot_ini, 0*Pf_react2, 'o', 'MarkerSize', 7.5, ...
                'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
            plot(rhot_lb, -(Pf_lb-Pf_react2)/(Pf_lb-Pf_react2), 'o', 'MarkerSize', 7.5, ...
                'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0]);

            ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$', 'FontSize', FS)
            xlabel('$\rho_\mathrm{total}/\rho_\mathrm{fluid}$', 'FontSize', FS);

            if icase == 1
                title('\textbf{(d) Thermodynamic relation}', ...
                    'Interpreter', 'latex', 'FontSize', FS);
            elseif icase==2
                title('\textbf{(e)}', 'FontSize', FS);
            elseif icase==3
                title('\textbf{(f)}', 'FontSize', FS);
            end

            if icase == 1
                text(0.50, 0.40, '\textbf{Hydraulic diffusion}', ...
                    'Units','normalized', 'Interpreter','latex', ...
                    'FontSize', FS, 'Color', col_f0_num);
            elseif icase == 2
                text(0.65, 0.8, '\textbf{1 reaction}', ...
                    'Units','normalized', 'Interpreter','latex', ...
                    'FontSize', FS, 'Color', col_f01_num);
            elseif icase == 3
                text(0.65, 0.8, '\textbf{1st reaction}', ...
                    'Units','normalized', 'Interpreter','latex', ...
                    'FontSize', FS, 'Color', col_f1_num);
                text(0.2, 0.6, '\textbf{2nd reaction}', ...
                    'Units','normalized', 'Interpreter','latex', ...
                    'FontSize', FS, 'Color', col_f2_num);
            end

            if icase == 3
                legend('Initial condition', 'Boundary condition', ...
                    'Location','southeast', ...
                    'Interpreter','latex', ...
                    'FontSize', FS);
            end

            axis([0.3 4.3 -1.1 0.1]);
            grid on;
            set(gca,'FontSize',FS);

            % Pressure profile
            subplot(3,2,1 + (icase-1)*2); box on; hold off;
            plot(x, -(Pf0-Pf_react2)/(Pf_lb-Pf_react2), 'k--', x, -(Pf-Pf_react2)/(Pf_lb-Pf_react2), 'k-', 'LineWidth', LW_PROFILE, 'HandleVisibility', 'off');

            if icase == 1
                xline(xf1_ana(end), '-',  'Color', col_f0_num,  'LineWidth', LW_FRONT);
            elseif icase == 2
                xline(xf1_ana_last,  '--', 'Color', col_f0_num,  'LineWidth', LW_FRONT);
                xline(xf01_num(end), '-',  'Color', col_f01_num, 'LineWidth', LW_FRONT);
            elseif icase == 3
                xline(xf1_ana_last,  '--', 'Color', col_f0_num,  'LineWidth', LW_FRONT);
                xline(xf01_num_last, '--', 'Color', col_f01_num, 'LineWidth', LW_FRONT);
                xline(xf2_num(end),  '-',  'Color', col_f1_num,  'LineWidth', LW_FRONT);
                xline(xf1_num(end),  '-',  'Color', col_f2_num,  'LineWidth', LW_FRONT);
            end

            ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$', 'FontSize', FS)
            xlabel('$x/L_c$', 'FontSize', FS);
            if icase == 1
            title('\textbf{(a) Pressure profile}', 'FontSize', FS);
            elseif icase==2
            title('\textbf{(b)}', 'FontSize', FS);
            elseif icase==3
            title('\textbf{(c)}', 'FontSize', FS);
            end

            lgd = legend( ...
                'Hydraulic diffusion', ...
                '1 reaction', ...
                '2 reactions, 1st front', ...
                '2 reactions, 2d front', ...
                'Location','east', ...
                'Interpreter','latex', ...
                'FontSize', FS-2);

            title(lgd, ['\textbf{Front position}: $t/t_c = ', num2str(t_pert,3), '$'], ...
                'Interpreter', 'latex');

            grid on;
            set(gca,'FontSize',FS);

            hold on;
            plot(0, -(Pf_lb-Pf_react2)/(Pf_lb-Pf_react2), 'o', 'MarkerSize', 7.5, ...
                'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0], ...
                'HandleVisibility', 'off');
            plot(x(end), -(Pf_react2-Pf_react2)/(Pf_lb-Pf_react2), 'o', 'MarkerSize', 7.5, ...
                'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
                'HandleVisibility', 'off');

            drawnow
        end
    end
    pause(2)
end

annotation('textbox',[0.195 0.94 0.64 0.03], ...
    'String','\textbf{Hydraulic diffusion}', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',FS*1.3);

annotation('textbox',[0.195 0.63 0.64 0.03], ...
    'String','\textbf{Hydraulic diffusion + 1 reaction}', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',FS*1.3);

annotation('textbox',[0.195 0.330 0.64 0.03], ...
    'String','\textbf{Hydraulic diffusion + 2 reactions}', ...
    'Interpreter','latex', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',FS*1.3);

saveas(gcf, fullfile(fig_dir, 'FIG_6.png'));
