function [m,s] = plotFy(daqwin)

nv = [142 143 144 145];
rv = [13 19 30 37];

for i = 1:4
    
    n = nv(i);
    r = rv(i);
    
%     figure(2)
%     plot(daqwin(:,174),daqwin(:,n),'b')
%     hold on
%     
%     figure(2)
%     plot(daqwin(find(daqwin(:,n)<9999),174),daqwin(find(daqwin(:,n)<9999),n),'--r')
%     hold on
%     
%     figure(200)
%     errorbar(r,mean(daqwin(find(daqwin(:,n)<9999),n)),std(daqwin(find(daqwin(:,n)<9999),n)),'ob')
%     hold on
    
    m(i) = mean(daqwin(find(daqwin(:,n)<9999),n));
    s(i) = std(daqwin(find(daqwin(:,n)<9999),n));
    
end

hFig = figure(1);
set(hFig, 'Position', [100 100 1200 800])

figure(1)
subplot(2,2,1)
errorbar(rv,m,s) 
ylabel('Fy')
% plot(rv,m,'-dk','DisplayName','exp.')
hold on
grid on
legend('-DynamicLegend','Location','southeast');

figure(1)
subplot(2,2,2)
% plot(rv,s,'-dk')
hold on
grid on

hFig = figure(10);
set(hFig, 'Position', [100 100 1200 800])

figure(10)
subplot(2,2,1)
plot(rv,m,'-dk','DisplayName','exp.')
ylabel('Fy - mean')
hold on
grid on
legend('-DynamicLegend','Location','southeast');

figure(10)
subplot(2,2,2)
plot(rv,s,'-dk')
ylabel('Fy - std')
hold on
grid on

end