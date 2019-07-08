%% plot MC samples:
figure


plot(X_ED2(:,1),X_ED2(:,2),'x')
xlim([-10 10])
ylim([-10 10])
grid
xlabel('q_{1}');
ylabel('q_{2}');

set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)

fig_name = 'linear_portfolio_q1q2';

export_fig(fig_name,'-png','-transparent');

%% plot output as function of x1:
figure
plot(X_ED2(:,1),Y_ED2,'x')

xlim([-10 10])
ylim([-10 10])
grid
xlabel('q_{1}');
ylabel('Y');

set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)


fig_name = 'linear_portfolio_q1y';

export_fig(fig_name,'-png','-transparent');

%% plot output as function of x2:
figure
plot(X_ED2(:,2),Y_ED2,'x')

xlim([-10 10])
ylim([-10 10])
grid
xlabel('q_{2}');
ylabel('Y');

set(gcf,'Color','w')
set(gca,'LineWidth',1)
set(gca,'FontSize',14)


fig_name = 'linear_portfolio_q2y';

export_fig(fig_name,'-png','-transparent');