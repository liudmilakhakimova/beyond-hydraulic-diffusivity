% FIG_3  Numerical and analytical propagation of two reaction fronts.
%
% Solves one-dimensional total mass conservation with Darcy flux using
% explicit time stepping. Pressure is obtained from a piecewise-linear
% pressure-density relation with two reaction intervals. Devolatilization
% and volatilization are calculated in separate cases.
%
% Compares numerical front positions with the coupled quasi-stationary
% solution in Appendix A. Panels show the constitutive relation, pressure,
% total density and front-position histories. Local front 1 is trailing
% and local front 2 is leading, opposite to the numbering in Appendix A.
%
% Input: physical and numerical parameters specified in this script.
% Output: figures/FIG_3.png and figures/DEVOL_2Fronts.gif.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear,figure(1),clf,colormap jet

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

file_name   = fullfile(fig_dir, 'DEVOL_2Fronts.gif');
irec        = 0;

for icase = 1:2
    % Physical parameters
    k_etaf        = 1;          % permeability/fluid viscosity, m2/Pa/s
    Lx            = 1;          % model length in x direction, m
    rhof          = 1;          % fluid density, kg/m^3
    t_pert        = 2.5;       % time of the experiment, s
    % Constitutive relation
    rhotp         = [   2.3   4     ]; % total density after the 1st and 2nd reaction
    rhotm         = [   1     2.4   ]; % total density before the 1st and 2nd reaction
    beta          = [5e1    1     1 ]; % compressibility in the 1st, 2d and 3d assemblages, 1/Pa
    Pf_react1     = 1.0;                                                  % pressure of the first reaction
    Pf_react2     = Pf_react1 + (rhotm(2) - rhotp(1))/(rhotm(2)*beta(2)); % pressure of the second reaction
    % Piecewise-linear pressure-density relation
    rhot_min = 0.1*rhotm(1);
    rhot_max = 1.05*rhotp(2);
    lutN     = 40000;
    rhot_lut = linspace(rhot_min,rhot_max,lutN);
    Pf_lut   = zeros(size(rhot_lut));
    dPdrho = [1/(rhotm(1)*beta(1)) 1/(rhotm(2)*beta(2)) 1/(rhotp(2)*beta(3))];
    for i = 1:length(rhot_lut)
        rho = rhot_lut(i);
        if rho <= rhotm(1)
            % first assemblage
            Pf_lut(i) = Pf_react1 + (rho-rhotm(1))*dPdrho(1);
        elseif rho <= rhotp(1)
            % first reaction
            Pf_lut(i) = Pf_react1;
        elseif rho <= rhotm(2)
            % second assemblage
            Pf_lut(i) = Pf_react1 + (rho-rhotp(1))*dPdrho(2);
        elseif rho <= rhotp(2)
            % second reaction
            Pf_lut(i) = Pf_react2;
        else
            % third assemblage
            Pf_lut(i) = Pf_react2 + (rho-rhotp(2))*dPdrho(3);
        end
    end
    % Numerical parameters
    nx   = 250;
    nt   = 3e5;
    nout = 500;
    dx   = Lx/(nx-1);
    x    = linspace(0,Lx,nx);
    dt   = dx^2/(rhof*k_etaf*max(dPdrho))/2.1;
    % Initial and boundary conditions
    if icase == 1
        rhot_lb   = 0.3*rhotm(1);    % total density at the left boundary
        rhot_ini  = 1.00*rhotp(2);   % total density of the initial assemblage
        Pf_lb     = interp1(rhot_lut,Pf_lut,rhot_lb ,'linear','extrap'); % pressure at the left boundary
        rhot      = rhot_ini*ones(size(x));
        rhot(1)   = rhot_lb;
    else
        rhot_lb   = 1.03*rhotp(2);    % total density at the left boundary
        rhot_ini  = 1.00*rhotm(1);   % total density of the initial assemblage
        Pf_lb     = interp1(rhot_lut,Pf_lut,rhot_lb ,'linear','extrap'); % pressure at the left boundary
        rhot      = rhot_ini*ones(size(x));
        rhot(1)   = rhot_lb;
    end
    % Local front 1 is trailing; local front 2 is leading.
    % Appendix A uses the opposite front numbering.
    % Quasi-stationary two-front solution
    if icase == 1
        DeltaRhot = [rhotp(1) - rhotm(1), rhotp(2) - rhotm(2)]; % total density jump due to 1st and 2nd reaction
        DeltaPf   = [Pf_react1-Pf_lb, Pf_react2-Pf_react1];     % pressure difference before 1st and 2nd reaction front
        Liu1  = DeltaRhot(1)/(rhotm(1)*beta(1)*DeltaPf(1));
        Liu2  = DeltaRhot(2)/(rhotm(2)*beta(2)*DeltaPf(2));
        Deff1 = rhof*k_etaf/(rhotm(1)*beta(1)*Liu1);            % Effective diffusivity associated with the 1st reaction front
        Deff2 = rhof*k_etaf/(rhotm(2)*beta(2)*Liu2);            % Effective diffusivity associated with the 2nd reaction front
    else
        DeltaRhot = [rhotp(2) - rhotm(2), rhotp(1) - rhotm(1)]; % total density jump due to 1st and 2nd reaction
        DeltaPf   = [abs(Pf_react2-Pf_lb), abs(Pf_react2-Pf_react1)];     % pressure difference before 1st and 2nd reaction front
        Liu1  = DeltaRhot(1)/(rhotp(2)*beta(3)*DeltaPf(1));
        Liu2  = DeltaRhot(2)/(rhotp(1)*beta(2)*DeltaPf(2));
        Deff1 = rhof*k_etaf/(rhotp(2)*beta(3)*Liu1);            % Effective diffusivity associated with the 1st reaction front
        Deff2 = rhof*k_etaf/(rhotp(1)*beta(2)*Liu2);            % Effective diffusivity associated with the 2nd reaction front
    end
    beta_qs = DeltaRhot(2)/DeltaRhot(1);
    r = (-(Deff1 + beta_qs*Deff2) + sqrt((Deff1 + beta_qs*Deff2)^2 + 4*Deff1*Deff2))/(2*Deff2);
    lambda2 = sqrt(2*Deff2/(1-r));
    lambda1 = r*lambda2;
    time    = [];
    xf1_num = [];
    xf2_num = [];
    xf1_ana = [];
    xf2_ana = [];
    it      = 0;
    while it*dt <= t_pert
        it = it + 1;
        % pressure from LUT
        Pf = interp1(rhot_lut,Pf_lut,rhot,'linear','extrap');
        % Darcy flux + total mass conservation
        qD = -k_etaf*diff(Pf)/dx;
        rhot(2:end-1) = rhot(2:end-1) - dt*diff(rhof*qD)/dx;
        % Reaction-front positions
        if mod(it,nout)==0 || it==1

            irec        = irec + 1;

            % numerical
            time(end+1) = it*dt;
            if icase == 1
                idx1 = find(Pf <= Pf_react1,1,'last');
                xf1_num(end+1) = interp1(Pf(idx1:idx1+1),x(idx1:idx1+1),Pf_react1,'linear','extrap');
                idx2 = find(Pf <= Pf_react2,1,'last');
                xf2_num(end+1) = interp1(Pf(idx2:idx2+1),x(idx2:idx2+1),Pf_react2,'linear','extrap');
            else
                idx1 = find(Pf >= Pf_react2,1,'last');
                xf1_num(end+1) = interp1(Pf(idx1:idx1+1),x(idx1:idx1+1),Pf_react2,'linear','extrap');
                idx2 = find(Pf >= Pf_react1,1,'last');
                xf2_num(end+1) = interp1(Pf(idx2:idx2+1),x(idx2:idx2+1),Pf_react1,'linear','extrap');
            end
            % analytical
            xf1_ana(end+1) = lambda1*sqrt(time(end));
            xf2_ana(end+1) = lambda2*sqrt(time(end));
            col_f1_num = [0.00 0.35 0.70];   % dark blue
            col_f2_num = [0.75 0.15 0.10];   % dark red
            col_f1_ana = [0.30 0.60 0.90];   % light blue
            col_f2_ana = [0.90 0.45 0.35];   % light red
            if icase == 1
                figure(1); clf;
            else
                figure(1);
            end
            set(gcf,'Color','w');

            if icase==1;
                subplot('position',[0.07 0.78 0.36 0.16]);
            else
                subplot('position',[0.5 0.78 0.36 0.16]);
            end
            box on;
            if icase == 1
                plot(rhot_lut, (Pf_lut-1.04167)/abs(min(Pf-1.04167)), 'k-', 'LineWidth', 2.0); hold on
                plot(rhot, (Pf-1.04167)/abs(min(Pf-1.04167)), 'bo', 'MarkerSize', 7.5, 'LineWidth', 1.2);
                plot(rhot_ini, (Pf_react2-1.04167)/abs(min(Pf-1.04167)), 'ko', 'MarkerSize', 7.5, 'LineWidth', 1.2);
                plot(rhot_lb, (Pf_lb-1.04167)/abs(min(Pf-1.04167)), 'ro', 'MarkerSize', 7.5, 'LineWidth', 1.2);
            else
                plot(rhot_lut, (Pf_lut-min(Pf))/abs(max(Pf-min(Pf))), 'k-', 'LineWidth', 2.0); hold on
                plot(rhot, (Pf-min(Pf))/abs(max(Pf-min(Pf))), 'bo', 'MarkerSize', 7.5, 'LineWidth', 1.2);
                plot(rhot_ini, (Pf_react1-min(Pf))/abs(max(Pf-min(Pf))), 'ko', 'MarkerSize', 7.5, 'LineWidth', 1.2);
                plot(rhot_lb, (Pf_lb-min(Pf))/abs(max(Pf-min(Pf))), 'ro', 'MarkerSize', 7.5, 'LineWidth', 1.2);
            end
            ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$', 'FontSize', 14)
            xlabel('$\rho_\mathrm{total}/\rho_\mathrm{fluid}$', 'FontSize', 14);
            if icase == 1
                title(['\textbf{(a) Devolatilization:} $t/t_c = ', num2str(it*dt/t_pert,3), '$'], 'FontSize', 14);
            else
                title(['\textbf{(e) Volatilization:} $t/t_c = ', num2str(it*dt/t_pert,3), '$'], 'FontSize', 14);
            end
            legend( ...
                'Thermodynamic relation', ...
                'Values inside model', ...
                'Initial value', ...
                'Boundary value', ...
                'Location','northwest', ...
                'Interpreter','latex', ...
                'FontSize',12.5);
            grid on;
            hold off

            if icase==1;subplot('position',[0.07 0.54 0.36 0.16]);else subplot('position',[0.5 0.54 0.36 0.16]); end
            box on;
            if icase==1
                plot(x, (Pf-1.04167)/abs(min(Pf-1.04167)), 'k-', 'LineWidth', 2.0);
            else
                plot(x, (Pf-min(Pf))/abs(max(Pf-min(Pf))), 'k-', 'LineWidth', 2.0);
            end
            xline(xf1_num(end), '--', 'Color', col_f1_num, 'LineWidth', 1.2);
            xline(xf2_num(end), '--', 'Color', col_f2_num, 'LineWidth', 1.2);
            ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$', 'FontSize', 14)
            xlabel('$x/L_c$', 'FontSize', 14);
            title('\textbf{(b) Pressure}', 'FontSize', 14);
            if icase == 2; title('\textbf{(f) Pressure}', 'FontSize', 14); end
            legend( ...
                '$P$', ...
                '$x_{\mathrm{front},2}$', ...
                '$x_{\mathrm{front},1}$', ...
                'Location','east', ...
                'Interpreter','latex', ...
                'FontSize',14);
            grid on;

            if icase==1;subplot('position',[0.07 0.30 0.36 0.16]);else subplot('position',[0.5 0.30 0.36 0.16]); end
            box on;
            plot(x, rhot, 'k-', 'LineWidth', 2.0);
            xline(xf1_num(end), '--', 'Color', col_f1_num, 'LineWidth', 1.2);
            xline(xf2_num(end), '--', 'Color', col_f2_num, 'LineWidth', 1.2);
            ylabel('$\rho_\mathrm{total}/\rho_\mathrm{fluid}$', 'FontSize', 14);
            xlabel('$x/L_c$', 'FontSize', 14);
            title('\textbf{(c) Total density}', 'FontSize', 14);
            if icase == 2; title('\textbf{(g) Total density}', 'FontSize', 14); end
            legend( ...
                '$\rho_{\mathrm{total}}$', ...
                '$x_{\mathrm{front},2}$', ...
                '$x_{\mathrm{front},1}$', ...
                'Location','east', ...
                'Interpreter','latex', ...
                'FontSize',14);
            grid on;

            if icase==1;subplot('position',[0.07 0.06 0.36 0.16]);else subplot('position',[0.5 0.06 0.36 0.16]); end
            hold on; box on;
            plot(time/t_pert, xf1_num, '-',  'Color', col_f1_num, 'LineWidth', 2.2);
            plot(time/t_pert, xf2_num, '-',  'Color', col_f2_num, 'LineWidth', 2.2);
            plot(time/t_pert, xf1_ana, '--', 'Color', col_f1_ana, 'LineWidth', 2.2);
            plot(time/t_pert, xf2_ana, '--', 'Color', col_f2_ana, 'LineWidth', 2.2);
            hold off
            xlabel('$t/t_c$', 'FontSize', 14);
            ylabel('$x_{\mathrm{front}}/L_c$', 'FontSize', 14);
            if icase == 1
                title('\textbf{(d) Devolatilization front}', 'FontSize', 14);
            else
                title('\textbf{(h) Volatilization front}', 'FontSize', 14);
            end
            legend({ ...
                '$x_{\mathrm{front},2}$ numerical', ...
                '$x_{\mathrm{front},1}$ numerical', ...
                '$x_{\mathrm{front},2}$ analytical', ...
                '$x_{\mathrm{front},1}$ analytical'}, ...
                'Location','northwest', ...
                'Interpreter','latex', ...
                'FontSize',14);
            grid on;

            set(gcf,"Position",[101 78 1025 788])
            drawnow

            [imind,cm]  = rgb2ind(frame2im(getframe(1)),128);
            if     irec == 1
                imwrite(imind,cm,file_name,'gif','Loopcount',inf,'DelayTime',2.0);
            else
                imwrite(imind,cm,file_name,'gif','WriteMode','append','DelayTime',0.05);
            end

        end
    end
end

imwrite(imind,cm,file_name,'gif','WriteMode','append','DelayTime',3.0);

saveas(gcf, fullfile(fig_dir, 'FIG_3.png'));
