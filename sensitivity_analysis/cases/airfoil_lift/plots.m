%% plot response surface

xlabel('C_{L} [-]');
ylabel('V [m/s]');
zlabel('L [N]');

set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)

fig_name = 'airfoil_lift_response';

export_fig(fig_name,'-png','-transparent');



