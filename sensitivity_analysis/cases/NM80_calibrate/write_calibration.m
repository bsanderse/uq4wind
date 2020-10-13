%% Write calibrated files
% This file reads the refernce polars and writes the calibrated polars.
% For the sake of completeness, the calibrated polars are plotted alongside.
% Note: The calibrated files for including the headerlines (necessary for Aero-Module run)
% is written in cal_files.m

addpath('../../AEROmodule/NM80_calibrate')

sec03 = readtable('section03_ref.dat');
p1 = table(sec03.Var1, sec03.Var2*(1+Delta_1), sec03.Var3, sec03.Var4);
writetable(p1, '3c.dat', 'Delimiter',' ', 'WriteVariableNames',0);

sec05 = readtable('section05_ref.dat');
p2 = table(sec05.Var1, sec05.Var2*(1+Delta_2), sec05.Var3, sec05.Var4);
writetable(p2, '5c.dat', 'Delimiter',' ', 'WriteVariableNames',0);

sec08 = readtable('section08_ref.dat');
p3 = table(sec08.Var1, sec08.Var2*(1+Delta_3), sec08.Var3, sec08.Var4);
writetable(p3, '8c.dat', 'Delimiter',' ', 'WriteVariableNames',0);

sec10 = readtable('section10_ref.dat');
p4 = table(sec10.Var1, sec10.Var2*(1+Delta_4), sec10.Var3, sec10.Var4);
writetable(p4, '10c.dat', 'Delimiter',' ', 'WriteVariableNames',0);

%% Plotting

% Fy03
figure()
plot (sec03.Var1, sec03.Var2,'k-','LineWidth',1.5)
hold on
plot (p1.Var1, p1.Var2,'r-','LineWidth',1.5)
xlabel('AOA[^{\circ}]')
ylabel('C_{L} [-]')
xlim([0 50])
title('Fy03')
legend('Reference', 'Calibrated','Location','southeast')
grid on

% Fy05
figure()
plot (sec05.Var1, sec05.Var2,'k-','LineWidth',1.5)
hold on
plot (p2.Var1, p2.Var2,'r-','LineWidth',1.5)
xlabel('AOA[^{\circ}]')
ylabel('C_{L} [-]')
xlim([0 50])
title('Fy05')
legend('Reference', 'Calibrated','Location','southeast')
grid on

% Fy08
figure()
plot (sec08.Var1, sec08.Var2,'k-','LineWidth',1.5)
hold on
plot (p3.Var1, p3.Var2,'r-','LineWidth',1.5)
xlabel('AOA[^{\circ}]')
ylabel('C_{L} [-]')
xlim([0 50])
title('Fy08')
legend('Reference', 'Calibrated','Location','southeast')
grid on


% Fy10
figure()
plot (sec10.Var1, sec10.Var2,'k-','LineWidth',1.5)
hold on
plot (p4.Var1, p4.Var2,'r-','LineWidth',1.5)
xlabel('AOA[^{\circ}]')
ylabel('C_{L} [-]')
xlim([0 50])
title('Fy10')
legend('Reference', 'Calibrated','Location','southeast')
grid on
hold off

