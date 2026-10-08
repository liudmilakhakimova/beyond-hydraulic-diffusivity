% FIG_4  Velocity ratio of two coupled reaction fronts.
%
% Evaluates the coupled-front solution over a range of effective-diffusivity
% and total-density-jump ratios. The positive root in equation (A31) gives
% r = v2/v1 = lambda2/lambda1, with front 1 leading and front 2 trailing.
%
% Deff1 is set to unity, so Deff2_vec is the diffusivity ratio Deff2/Deff1.
% gam is Delta_rho_total,1/Delta_rho_total,2; the vertical axis shows 1/gam.
% The red marker denotes the representative gypsum-experiment parameters
% discussed in section 4.4 of the manuscript.
%
% Input: parameter ranges and marker coordinates specified in this script.
% Output: figures/FIG_4.png.
%
% Companion code for:
% Beyond hydraulic diffusion: reaction front propagation and
% timescales in metamorphic (de)volatilization processes.
% Liudmila Khakimova, Stefan M. Schmalholz and Yury Y. Podladchikov.
% University of Lausanne. Contact: liudmila.khakimova@unil.ch

clear all, close all, clc, colormap(jet)

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% r = v2/v1, with front 1 leading and front 2 trailing.
Deff1           = 1;
Deff2_vec       = 10.^[-8:0.1:2];
gam             = 10.^[-2:0.1:2];

% Coupled-front solution (Appendix A).
for it=1:length(Deff2_vec)
    Deff2       = Deff2_vec(it);    % effective diffusivity for water-vapor reaction front
    r           = (-(Deff2+gam.*Deff1) + sqrt((Deff2+gam.*Deff1).^2+4.*Deff1.*Deff2))./(2*Deff1);
    rMat(it,:)  = r;
    Lam1        = sqrt(2.*Deff1./(1-r));
    Lam2        = r.*Lam1;
    v2_v1(it,:) = Lam2./Lam1;
end

subplot(221)
contourf(Deff2_vec,1./gam,(v2_v1)',[1e-7 1e-6 1e-5 1e-4 0.001 0.01:0.01:0.09 0.1:0.05:1],'edgecolor','none'), hold on
[c,h] = contour(Deff2_vec,1./gam,(v2_v1)',[0.4 0.6 0.8 0.9 0.98],'-w','LineWidth',2);
clabel(c,h,'color','k','fontsize',14,'LabelSpacing',195)
[c,h] = contour(Deff2_vec,1./gam,(v2_v1)',[0.0001 0.001 0.01 0.1 0.2],'-w','LineWidth',2);
clabel(c,h,'color','w','fontsize',14,'LabelSpacing',320)

plot(5e-5,3.2,'+r','MarkerSize',8,'LineWidth',2)

set(gca,'xscale','log','yscale','log')
caxis([0 1])
colorbar
set(gca,'FontSize',16)
ylabel('$\Delta \rho_\mathrm{total,2} / \Delta \rho_\mathrm{total,1}$','FontSize',18)
xlabel('$D_\mathrm{{eff,2}} / D_\mathrm{{eff,1}}$','FontSize',18)
title('Contours of $v_2 / v_1$','FontSize',18)
axis([4e-5 1e2 1e-1 10])

set(gcf,'Position', [210 73 1133 793])
saveas(gcf, fullfile(fig_dir, 'FIG_4.png'));
