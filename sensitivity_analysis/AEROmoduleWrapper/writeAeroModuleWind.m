function writeAeroModuleWind(u,folder)

filename = [folder,'\wind.dat'];
fid = fopen(filename,'w');
fprintf(fid,'!time [s]  u [m/s]  v [m/s]  w [m/s]\n');
fprintf(fid,'  0.0      %f    0.0      0.0\n',u);
fprintf(fid,'  5.0      %f    0.0      0.0\n',u);
fprintf(fid,' 10.0      %f    0.0      0.0\n',u);
fprintf(fid,' 20.0      %f    0.0      0.0\n',u);
fclose(fid);

end