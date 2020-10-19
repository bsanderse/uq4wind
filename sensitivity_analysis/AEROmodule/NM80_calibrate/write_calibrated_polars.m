%% Adding headerlines to calibrated polars

% Fy03
fidin1 = fopen('3c.dat', 'rt');
fidout1 = fopen('section03_ref_c.dat', 'wt');
fprintf(fidout1, '%s\n', '! Aero module input file for airfoil data',...
    '', 'Airfoil_Name  Section03', 't/c            0.333000         ! thickness ratio w.r.t. chord',...
    '','format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm','',...
    'Reynolds_Nr 10000000.000000');
while true
  thisline = fgets(fidin1);
  if ~ischar(thisline); break; end   %end of file
  fwrite(fidout1, thisline);
end
fclose(fidout1);
fclose(fidin1);

% Fy05

fidin2 = fopen('5c.dat', 'rt');
fidout2 = fopen('section05_ref_c.dat', 'wt');
fprintf(fidout2, '%s\n', '! Aero module input file for airfoil data',...
    '', 'Airfoil_Name  Section05', 't/c            0.243000         ! thickness ratio w.r.t. chord',...
    '','format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm','',...
    'Reynolds_Nr 10000000.000000');
while true
  thisline = fgets(fidin2);
  if ~ischar(thisline); break; end   %end of file
  fwrite(fidout2, thisline);
end
fclose(fidout2);
fclose(fidin2);

% Fy08

fidin3 = fopen('8c.dat', 'rt');
fidout3 = fopen('section08_ref_c.dat', 'wt');
fprintf(fidout3, '%s\n', '! Aero module input file for airfoil data',...
    '', 'Airfoil_Name  Section08', 't/c            0.197000         ! thickness ratio w.r.t. chord',...
    '','format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm','',...
    'Reynolds_Nr 10000000.000000');
while true
  thisline = fgets(fidin3);
  if ~ischar(thisline); break; end   %end of file
  fwrite(fidout3, thisline);
end
fclose(fidout3);
fclose(fidin3);

% Fy10

fidin4 = fopen('10c.dat', 'rt');
fidout4 = fopen('section10_ref_c.dat', 'wt');
fprintf(fidout4, '%s\n', '! Aero module input file for airfoil data',...
    '', 'Airfoil_Name  Section10', 't/c            0.187000         ! thickness ratio w.r.t. chord',...
    '','format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm','',...
    'Reynolds_Nr 10000000.000000');
while true
  thisline = fgets(fidin4);
  if ~ischar(thisline); break; end   %end of file
  fwrite(fidout4, thisline);
end
fclose(fidout4);
fclose(fidin4);


