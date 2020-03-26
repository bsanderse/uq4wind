function [daqwin] = filt(daqwin)

t = daqwin(:,174);

% figure(1)
% plot(daqwin(:,174),daqwin(:,84),'+b')
% hold on
% grid on

index = find(t>=200 & t<=450);

% figure(1)
% plot(daqwin(index,174),daqwin(index,84),'or')
% hold on
% grid on

daqwinfilt = daqwin(index,:);

daqwinfilt(:,174) = daqwinfilt(:,174)-daqwinfilt(1,174);

clear daqwin

daqwin = daqwinfilt;

end