%% Write calibrated files
% This file reads the reference polars and writes the calibrated polars.
% For the sake of completeness, the calibrated polars are plotted alongside.
% Note: The calibrated files for including the headerlines (necessary for Aero-Module run)
% is written in cal_files.m

addpath('AEROmodule/NM80_calibrate')

% repeat the steps from writeAeroModuleInput to get the perturbed lift
% curves
n_polar = P{31}{1};
aoa     = P{31}{6}; % angle of attack
ind_aoa = P{31}{7+n_polar}; % indices that are perturbed
d       = ones(length(ind_aoa),1);
ind_CL  = 1;
plotCurve = 1;

for i=1:4
    X_CL    = BayesianAnalysis.Results.PostProc.PointEstimate.X(i);
    CL_PERTURB = P{26}{i}{3};
    CL_REF     = P{31}{6+i}{ind_CL};
    CL{i}   = computeCurves(1, ind_aoa, X_CL*d, CL_PERTURB*d, plotCurve, ...
                aoa, CL_REF, 3, 1:length(aoa));    
end


%%
% sec03 = readtable('section03_ref.dat');
% p1 = table(sec03.Var1, sec03.Var2*(1+Delta(1)), sec03.Var3, sec03.Var4);
% writetable(p1, '3c.dat', 'Delimiter',' ', 'WriteVariableNames',0);
% 
% sec05 = readtable('section05_ref.dat');
% p2 = table(sec05.Var1, sec05.Var2*(1+Delta(2)), sec05.Var3, sec05.Var4);
% writetable(p2, '5c.dat', 'Delimiter',' ', 'WriteVariableNames',0);
% 
% sec08 = readtable('section08_ref.dat');
% p3 = table(sec08.Var1, sec08.Var2*(1+Delta(3)), sec08.Var3, sec08.Var4);
% writetable(p3, '8c.dat', 'Delimiter',' ', 'WriteVariableNames',0);
% 
% sec10 = readtable('section10_ref.dat');
% p4 = table(sec10.Var1, sec10.Var2*(1+Delta(4)), sec10.Var3, sec10.Var4);
% writetable(p4, '10c.dat', 'Delimiter',' ', 'WriteVariableNames',0);

%% Plotting

% Fy03
% figure()
% plot (sec03.Var1, sec03.Var2,'k-','LineWidth',1.5)
% hold on
% plot (p1.Var1, p1.Var2,'r-','LineWidth',1.5)
% xlabel('AOA[^{\circ}]')
% ylabel('C_{L} [-]')
% xlim([0 50])
% title('Fy03')
% legend('Reference', 'Calibrated','Location','southeast')
% grid on
% 
% % Fy05
% figure()
% plot (sec05.Var1, sec05.Var2,'k-','LineWidth',1.5)
% hold on
% plot (p2.Var1, p2.Var2,'r-','LineWidth',1.5)
% xlabel('AOA[^{\circ}]')
% ylabel('C_{L} [-]')
% xlim([0 50])
% title('Fy05')
% legend('Reference', 'Calibrated','Location','southeast')
% grid on
% 
% % Fy08
% figure()
% plot (sec08.Var1, sec08.Var2,'k-','LineWidth',1.5)
% hold on
% plot (p3.Var1, p3.Var2,'r-','LineWidth',1.5)
% xlabel('AOA[^{\circ}]')
% ylabel('C_{L} [-]')
% xlim([0 50])
% title('Fy08')
% legend('Reference', 'Calibrated','Location','southeast')
% grid on
% 
% 
% % Fy10
% figure()
% plot (sec10.Var1, sec10.Var2,'k-','LineWidth',1.5)
% hold on
% plot (p4.Var1, p4.Var2,'r-','LineWidth',1.5)
% xlabel('AOA[^{\circ}]')
% ylabel('C_{L} [-]')
% xlim([0 50])
% title('Fy10')
% legend('Reference', 'Calibrated','Location','southeast')
% grid on
% hold off

