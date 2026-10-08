% FIG_5  Ratio of effective to hydraulic diffusivity.
%
% Plots Deff/Dhyd = (rho_fluid/Delta_rho_total)*Delta_P/K, equation (25),
% for pressure differences from 0.1 to 100 MPa and K = 50 GPa.
%
% Vertical markers represent seven metamorphic reactions. Their density
% jumps are specified directly or calculated from solid densities, water
% mass fractions and conservation of non-volatile solid mass. The reaction
% parameters and literature references are listed below.
%
% Input: pressure range, modulus and reaction parameters in this script.
% Output: figures/FIG_5.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear all, close all, clc, figure(1)

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

subplot('position',[0.1 0.6 0.6 0.35])
rho_ratio   = [2:50];          % rho_fluid / Delta_rho_total
Delta_P     = [0.1 1 10 100]*1e6;   % Fluid pressure difference
K           = 50e9;             % Poroelastic modulus in hydraulic diffusion
col         = {'-c','-k','-r','b'};
for ip=1:length(Delta_P)
    Diff_Ratio  = rho_ratio.*Delta_P(ip)/K;
    hold on
    plot(rho_ratio,Diff_Ratio,col{ip},'LineWidth',2)
end
set(gca,'XScale','log','YScale','log','ytick',[1e-5 1e-4 1e-3 1e-2 1e-1])
set(gca,'FontSize',16)
title('Ratio of effective to hydraulic diffusivity')
xlabel('$\rho_\mathrm{fluid} / \Delta\rho_\mathrm{total}$','Interpreter','latex','FontSize',20)
ylabel('$D_\mathrm{eff} / D_\mathrm{hyd}$','Interpreter','latex','FontSize',20)
text(2.1,5e-2,'K = 50 GPa','Interpreter','latex','FontSize',16)
grid on

RR          = [950/(2700-2300) ];
PP          = [4.75e-3         ];

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

RR_L        = r_f ./ Delta_rho_tot_2;
PP_L        = 1e-4*ones(size(RR_L));
col         = {'-ok','-sk','-<k','-vk','-+k','-dk','-xk'};
for ida=1:length(RR_L)
    plot([RR_L(ida) RR_L(ida)],[3.5e-6 9e-2],col{ida},'LineWidth',0.5,'MarkerSize',10)
end
leg = legend('$\Delta$P = 0.1 MPa','$\Delta$P = 1 MPa','$\Delta$P = 10 MPa','$\Delta$P = 100 MPa',...
    '1) Dehy: gyp $\rightarrow$ bas', '2) Dehy: bru $\rightarrow$ per', '3) Dehy: atg+bru $\rightarrow$ ol',...
    '4) Dehy: serp $\rightarrow$ du', '5) Dehy: atg $\rightarrow$ fo+en',...
    '6) Hy: gra $\rightarrow$ ecl', '7) Hy: gra $\rightarrow$ ecl');
set(leg,'Interpreter','latex','FontSize',16,'Location','northeastoutside')

axis([2 50 3e-6 1e-1])
set(gca,'xtick',[2 4 8 10 20 30 40 50])

set(gcf,'position',[140 197 945 633])
saveas(gcf, fullfile(fig_dir, 'FIG_5.png'));
