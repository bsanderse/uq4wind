function [meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber)

figure(signalnumber)
plot(daqwin(:,174),daqwin(:,signalnumber),'b')
hold on
grid on

meanout = mean(daqwin(:,signalnumber));
stdout = std(daqwin(:,signalnumber));
minout = min(daqwin(:,signalnumber));
maxout = max(daqwin(:,signalnumber));

disp([name,': ',num2str(meanout),' ',num2str(stdout),' ',num2str(minout),' ',num2str(maxout)])

end