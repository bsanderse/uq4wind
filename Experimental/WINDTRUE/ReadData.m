function dat = ReadData(fp,nch)
fid = fopen(fp,'r');
dat=fread(fid,'single');
fclose(fid);
dat=reshape(dat,[nch length(dat)/nch])';
end