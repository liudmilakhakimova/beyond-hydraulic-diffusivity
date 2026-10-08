% FIG_2  Single-front dehydration and hydration calculations.
%
% Compares the saved numerical results of this study with those labelled
% Schmalholz et al. (2024). The left and right columns show dehydration
% and hydration, respectively. Each column contains the pressure-density
% relation, final pressure and density profiles, and front-position history.
%
% Black curves use DATA_DEHY_*.mat; red curves use DATA_Liudmilalike_*.mat.
% The blue dashed analytical histories are loaded from X_front_ana.
% Pressure, density, position and time retain their original normalization.
% The script plots saved results; it does not rerun the numerical solvers.
%
% Input: data/DATA_DEHY_Dehyd.mat, data/DATA_DEHY_Hyd.mat,
%        data/DATA_Liudmilalike_Dehyd.mat, data/DATA_Liudmilalike_Hyd.mat.
% Output: figures/FIG_2.png.
% Supplied comparison image: reference/FIG_2.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear variables, close all, clc

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

lw              = 1.5;
fs              = 14;

% Dehydration

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'DATA_DEHY_Dehyd.mat'));

subplot('position',[0.1 0.57 0.35 0.185])
plot(x,Pf,'-k','LineWidth',lw), hold on
subplot('position',[0.1 0.31 0.35 0.175])
plot(x,rhot,'-k','LineWidth',lw), hold on
subplot('position',[0.1 0.075 0.35 0.15])
plot(Time_vec,X_front,'-k','LineWidth',lw), hold on

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'DATA_Liudmilalike_Dehyd.mat'));

subplot('position',[0.1 0.85 0.35 0.115])
plot(rhot_lut,P_lut,'-k','LineWidth',lw), hold on
plot([1.913 ],[-1],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([1.9148],[ 0],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([     2],[ 0],'ko','MarkerSize',10,'LineWidth',lw*1.2)
ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$'), xlabel('$\rho_\mathrm{total} / \rho_\mathrm{fluid}$'), title('$\textbf{(a) Density vs fluid pressure}$')
set(gca,'fontsize',fs)
axis([1.9 2.02 -1 0.2])
grid on

subplot('position',[0.1 0.57 0.35 0.185])
plot(x,Pf,'-r','LineWidth',lw)
Pf_ini = 0*x; Pf_ini(1) = - dP;
plot(x,Pf_ini,'--k','LineWidth',lw)
plot([0],[-1],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([4.794],[ 0],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([4.854 6],[ 0 0],'ko','MarkerSize',10,'LineWidth',lw*1.2)
ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$'),xlabel('$x / L_c$'),title(['\textbf{(b) Pressure}; $t / t_c = ' num2str(it*dt),'$'])
legend('Numerics; Schmalholz et al. (2024)','Numerics; This study','Initial condition, both simulations','Location','southeast')
axis([-0.1 6 -1.1 0.1])
grid on
set(gca,'fontsize',fs)

subplot('position',[0.1 0.31 0.35 0.175])
plot(x,rhot,'-r','LineWidth',lw)
rhot_ini     = 0*x + rhotp(1);
rhot_ini(1)  = rhot_BC;
plot(x,rhot_ini,'--k','LineWidth',lw)
plot([0],[1.913],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([4.794],[ 1.9148],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([4.854 6],[ 2 2],'ko','MarkerSize',10,'LineWidth',lw*1.2)
plot([0.01 4.8 4.5 4.8 4.5]*0.99,[1.975 1.975 1.97 1.975 1.98],'LineWidth',lw/2,'color','k')
text(2,1.966,'$x_\mathrm{front}$','FontSize',14)
xlabel('$x / L_c$'),ylabel('$\rho_\mathrm{total} / \rho_\mathrm{fluid}$'),title(['$\textbf{(c) Total density}$'])
axis([-0.1 6 1.9 2.02])
grid on
set(gca,'fontsize',fs)

subplot('position',[0.1 0.075 0.35 0.15])
plot(Time_vec,X_front,'-r','LineWidth',lw)
plot(Time_vec(1:1e0:end),X_front_ana(1:1e0:end),'--b','LineWidth',lw)
legend('Numerics; Schmalholz et al. (2024)','Numerics; This study','Analytical solution','Location','southeast')
xlabel('$t / t_c$'), ylabel('$x_\mathrm{front} / L_c$'), title('$\textbf{(d) Dehydration front}$')
axis([0 1 0 5])
grid on
set(gca,'fontsize',fs)


% Hydration

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'DATA_DEHY_Hyd.mat'));

subplot('position',[0.535 0.57 0.35 0.185])
plot(x,Pf,'-k','LineWidth',lw), hold on
subplot('position',[0.535 0.31 0.35 0.175])
plot(x,rhot,'-k','LineWidth',lw), hold on
subplot('position',[0.535 0.075 0.35 0.15])
plot(Time_vec,X_front,'-k','LineWidth',lw), hold on

load(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'data', 'DATA_Liudmilalike_Hyd.mat'));

subplot('position',[0.535 0.85 0.35 0.115])
plot(rhot_lut,P_lut,'-k','LineWidth',lw), hold on
plot([2.2  ],[ 1],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([2.197],[ 0],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([2    ],[ 0],'ko','MarkerSize',10,'LineWidth',lw*1.2)
ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$'), xlabel('$\rho_\mathrm{total} / \rho_\mathrm{fluid}$'), title('$\textbf{(e) Density vs fluid pressure}$')
set(gca,'fontsize',fs)
axis([1.98 2.25 -0.2 1])
grid on

subplot('position',[0.535 0.57 0.35 0.185])
plot(x,Pf,'-r','LineWidth',lw)
Pf_ini = 0*x; Pf_ini(1) =  dP;
plot(x,Pf_ini,'--k','LineWidth',lw)
plot([0],[1],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([3.135],[ 0],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([3.196 6],[ 0 0],'ko','MarkerSize',10,'LineWidth',lw*1.2)
ylabel('$(P-P_\mathrm{amb}) / P_\mathrm{per}$'),xlabel('$x / L_c$'),title(['\textbf{(f) Pressure}; $t / t_c = ' num2str(it*dt),'$'])
legend('Numerics; Schmalholz et al. (2024)','Numerics; This study','Initial condition, both simulations')
axis([-0.1 6 -0.1 1.1])
grid on
set(gca,'fontsize',fs)
subplot('position',[0.535 0.31 0.35 0.175])
plot(x,rhot,'-r','LineWidth',lw)
rhot_ini     = 0*x + rhotp(1);
rhot_ini(1)  = rhot_BC;
plot(x,rhot_ini,'--k','LineWidth',lw)
plot([0],[2.2],'ro','MarkerSize',10,'LineWidth',lw*1.2)
plot([3.135],[ 2.197],'bo','MarkerSize',10,'LineWidth',lw*1.2)
plot([3.196 6],[ 2 2],'ko','MarkerSize',10,'LineWidth',lw*1.2)

xlabel('$x / L_c$'),ylabel('$\rho_\mathrm{total} / \rho_\mathrm{fluid}$'),title(['$\textbf{(g) Total density}$'])
axis([-0.1 6 1.95 2.25])
grid on
set(gca,'fontsize',fs)
subplot('position',[0.535 0.075 0.35 0.15])
plot(Time_vec,X_front,'-r','LineWidth',lw)
plot(Time_vec(1:1e0:end),X_front_ana(1:1e0:end),'--b','LineWidth',lw)
legend('Numerics; Schmalholz et al. (2024)','Numerics; This study','Analytical solution','Location','southeast')
xlabel('$t / t_c$'), ylabel('$x_\mathrm{front} / L_c$'), title('$\textbf{(h) Hydration front}$')
axis([0 1 0 3.5])
grid on
set(gca,'fontsize',fs)


set(gcf,'Position',[223 89 948 777])
saveas(gcf, fullfile(fig_dir, 'FIG_2.png'));
