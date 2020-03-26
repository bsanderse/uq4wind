function [m,s] = plotFx(daqwin)

nv = [138 139 140 141];
rv = [13 19 30 37];

for i = 1:4
    
    n = nv(i);
    r = rv(i);
    
%     figure(1)
%     plot(daqwin(:,174),daqwin(:,n),'b')
%     hold on
%     
%     figure(1)
%     plot(daqwin(find(daqwin(:,n)<9999),174),daqwin(find(daqwin(:,n)<9999),n),'--r')
%     hold on
%     
%     figure(100)
%     errorbar(r,mean(daqwin(find(daqwin(:,n)<9999),n)),std(daqwin(find(daqwin(:,n)<9999),n)),'ob')
%     hold on
    
    m(i) = mean(daqwin(find(daqwin(:,n)<9999),n));
    s(i) = std(daqwin(find(daqwin(:,n)<9999),n));
    
end

hFig = figure(1);
set(hFig, 'Position', [100 100 1200 800])

figure(1)
subplot(2,2,3)
errorbar(rv,m,s) 
ylabel('Fx')
hold on
grid on

figure(1)
subplot(2,2,4)
% plot(rv,s,'-dk')
hold on
grid on

hFig = figure(10);
set(hFig, 'Position', [100 100 1200 800])

figure(10)
subplot(2,2,3)
plot(rv,m,'-dk')
ylabel('Fx - mean')
hold on
grid on

figure(10)
subplot(2,2,4)
plot(rv,s,'-dk')
ylabel('Fx - std')
hold on
grid on

end