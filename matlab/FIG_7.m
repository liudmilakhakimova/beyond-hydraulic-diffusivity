% FIG_7  Permeability estimates and reaction-front propagation times.
%
% Panel (a) estimates permeability from the measured front coefficient
% x_front^2/t using equation (26), with prescribed viscosity, density ratio
% and pressure differences. Field examples are serpentinite alteration
% (Beinlich et al., 2020), blueschist metasomatism (John et al., 2012), and
% fluid-release veins (Taetz et al., 2018).
%
% Panel (b) calculates propagation time versus squared front distance for
% seven reaction-induced density changes, using the single-front solution
% in equation (17). The density calculation is the same as in FIG_5.
% Here x2 denotes squared distance, not the second front position.
%
% Input: field distances/times and material parameters in this script.
% Output: figures/FIG_7.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear all, close all, clc

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

subplot('position',[0.1 0.6 0.6 0.35])

x2_t        = 10.^[-11:0.1:-7];
rho_ratio   = 10;              % rho_fluid / Delta_rho_total
Delta_P     = [0.1 1 10 100]*1e6;   % Fluid pressure difference
eta         = 1e-4;             % Fluid viscosity
col         = {'-c','-k','-r','b'};

for ip=1:length(Delta_P)
    k           = eta/2/Delta_P(ip)/rho_ratio .* x2_t;
    hold on
    plot(x2_t,k,col{ip},'LineWidth',1.5)

end
set(gca,'XScale','log','YScale','log')
set(gca,'FontSize',16)
title('(a) Effective diffusivity versus permeability')
xlabel('$x_\mathrm{front}^2/t$ (m$^2$ s$^{-1}$)','Interpreter','latex','FontSize',20)
ylabel('$k$ (m$^2$)','Interpreter','latex','FontSize',20)

yr_2_s      = 365*24*3600;
% Ref         Fusseis   Beinlich    John2012    Taetz
xD          = [1e-3     2.6         1           0.04       ];
tD          = [1.2e4    22*yr_2_s   760*yr_2_s  0.18*yr_2_s];
x2_t        = xD.^2./tD;

col         = {'-ok','-sk','-dk','-xk'};
for ida=2:length(x2_t)
    plot([x2_t(ida) x2_t(ida)],[6e-25 9e-19],col{ida},'LineWidth',0.5,'MarkerSize',10)
end

axis([1e-11 2e-8 5e-25 1e-18])

leg = legend('$\Delta$P = 0.1 MPa','$\Delta$P = 1 MPa','$\Delta$P = 10 MPa','$\Delta$P = 100 MPa',...
    'Serpentinite alteration','Blueschist metasomatism','Fluid release veins');
set(leg,'Interpreter','latex','FontSize',16,'Location','northeastoutside')
text(3e-8,11e-24,'$\eta$ = 10$^{-4}$ Pa s','Interpreter','latex','FontSize',18)
text(3e-8,3e-24,'$\rho_{fluid} / \Delta\rho_{total}$ = 10','Interpreter','latex','FontSize',18)
grid on

subplot('position',[0.1 0.1 0.6 0.35])

% Reaction data; references follow the parameter arrays.
Rho_s_LP    = [3500  3200  3300  2685   3000  ];
Rho_s_HP    = [2275  2560  2492  3228   3200  ];
C_sH2O_LP   = [0     0     0     0      0.002 ];
C_sH2O_HP   = [0.31  0.17  0.126 0.007  0.007 ];
% (1) per - bru; Schmalholz et al. (2020)
% (2) oli - ser; Schmalholz et al. (2023)
% (3) Dunite - serp; Malvoisin et al., 2021
% (4) gra - ecl; Bras et al. (2023)
% (5) gra - ecl; Centrella (2019)
Rho_s_nophi     = [Rho_s_HP(1)  Rho_s_HP(2)  Rho_s_HP(3)  Rho_s_LP(4)  Rho_s_LP(5)  ];
Rho_s_phi       = [Rho_s_LP(1)  Rho_s_LP(2)  Rho_s_LP(3)  Rho_s_HP(4)  Rho_s_HP(5)  ];
C_nophi         = [C_sH2O_HP(1) C_sH2O_HP(2) C_sH2O_HP(2) C_sH2O_LP(4) C_sH2O_LP(5) ];
C_phi           = [C_sH2O_LP(1) C_sH2O_LP(2) C_sH2O_LP(2) C_sH2O_HP(4) C_sH2O_HP(5) ];
Phi_calc        = 1 - Rho_s_nophi.*(1-C_nophi) ./ (Rho_s_phi.*(1-C_phi));
r_f             = 1000;
Rho_tot_phi     = (1-Phi_calc).*Rho_s_phi   + Phi_calc*r_f;
Rho_tot_nophi   =               Rho_s_nophi;
Delta_rho_tot   = abs(Rho_tot_phi-Rho_tot_nophi);
Delta_rho_tot   = [88 Delta_rho_tot 52]; % add data for atg-out, negative Clapeyon slope; from Porkolab et al 2025
Delta_rho_tot_2 = [Delta_rho_tot(1:4) Delta_rho_tot(7) Delta_rho_tot(5:6)];

x2          = 10.^[-1:0.5:0];  % squared front distance, m^2
Delta_P     = 10*1e6;   % Fluid pressure difference
k           = [1e-21];
eta         = 1e-4;             % Fluid viscosity
col         = {'-sk','-xk','-ok','-+k','-dk','-sb','-ob'};
s_to_yr     = 1/(3600*24*365.25);
s_to_days   = 1/(3600*24);
s_to_months   = 1/(3600*24*12);
for ip=1:length(Delta_rho_tot_2)
    t           = eta/2/Delta_P/(r_f./Delta_rho_tot_2(ip)) .* (x2) /k;
    hold on
    plot(x2,t*s_to_yr,col{ip},'LineWidth',1.5,'MarkerSize',10)
end
set(gca,'XScale','log','YScale','log')
set(gca,'FontSize',16)
leg = legend('1) gyp$\rightarrow$bas; 88 kgm$^{-3}$', '2) bru$\rightarrow$per; 154 kgm$^{-3}$',...
    '3) atg+bru$\rightarrow$ol; 99 kgm$^{-3}$',...
    '4) serp$\rightarrow$du; 50 kgm$^{-3}$', '5) atg$\rightarrow$fo+en; 52 kgm$^{-3}$',...
    '6) gra$\rightarrow$ecl; 181 kgm$^{-3}$', '7) gra$\rightarrow$ecl; 73 kgm$^{-3}$');
set(leg,'Interpreter','latex','FontSize',16,'Location','northeastoutside')

text(1.2,2.0,'$\eta$ = 10$^{-4}$ Pa s','Interpreter','latex','FontSize',18)
text(1.2,1.4,'$\Delta P$ = 10 MPa','Interpreter','latex','FontSize',18)
text(1.2,1.0,'$k$ = 10$^{-21}$ m$^2$','Interpreter','latex','FontSize',18)
text(1.2,0.8,'$\rho_\mathrm{fluid}$ = 1000 kg m$^{-3}$','Interpreter','latex','FontSize',18)

title('(b) Time scale and density change')
xlabel('$x_\mathrm{front}^2$ (m$^2$)','Interpreter','latex','FontSize',20)
ylabel('$t$ (yr)','Interpreter','latex','FontSize',20)
grid on

set(gcf,'Position', [192 145 1071 616])
saveas(gcf, fullfile(fig_dir, 'FIG_7.png'));
