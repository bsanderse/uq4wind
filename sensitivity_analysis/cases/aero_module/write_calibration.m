%% Write calibrated files
sec03 = readtable(['C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section03_ref.dat']);
p1 = table(sec03.Var1, sec03.Var2*(1+Delta_1), sec03.Var3, sec03.Var4);
writetable(p1, 'C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section03_ref_c.dat', 'Delimiter',' ');

sec05 = readtable(['C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section05_ref.dat']);
p2 = table(sec05.Var1, sec05.Var2*(1+Delta_2), sec05.Var3, sec05.Var4);
writetable(p2, 'C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section05_ref_c.dat', 'Delimiter',' ');

sec08 = readtable(['C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section08_ref.dat']);
p3 = table(sec08.Var1, sec08.Var2*(1+Delta_3), sec08.Var3, sec08.Var4);
writetable(p3, 'C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section08_ref_c.dat', 'Delimiter',' ');

sec10 = readtable(['C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section10_ref.dat']);
p4 = table(sec10.Var1, sec10.Var2*(1+Delta_4), sec10.Var3, sec10.Var4);
writetable(p4, 'C:\Users\kriek\Desktop\windtrue\sensitivity_analysis\AEROmodule\NM80_calibrate/section10_ref_c.dat', 'Delimiter',' ');


