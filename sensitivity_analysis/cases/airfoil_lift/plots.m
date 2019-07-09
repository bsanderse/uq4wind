%% plot response surface

xlabel('C_{L} [-]');
ylabel('V [m/s]');
zlabel('L [N]');

%%
set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)

fig_name = 'airfoil_lift_inputs';

export_fig(fig_name,'-png','-transparent');

%% quadrature points


plot(myPCE_Quad.ExpDesign.X(:,1),myPCE_Quad.ExpDesign.X(:,2),'x')
xlabel('C_{L} [-]');
ylabel('V [m/s]');
set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)

fig_name = 'airfoil_lift_ED';

export_fig(fig_name,'-png','-transparent');


%% Sobol bars
fig_name = 'airfoil_lift_Sobol';

export_fig(fig_name,'-png','-transparent');

%% Sobol convergence

legend('Monte Carlo - CL','Monte Carlo - V','PCE - CL','PCE - V');

set(gcf,'Color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)

xlabel('number of simulations')

fig_name = 'airfoil_lift_Sobol_convergence';

export_fig(fig_name,'-png','-transparent');
