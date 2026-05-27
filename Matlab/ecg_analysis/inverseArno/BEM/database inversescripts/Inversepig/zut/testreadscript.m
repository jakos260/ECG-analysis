fname=mfilename;
fp=fopen([fname '.m'],'r');
codestr=fread(fp,inf,'char');
fclose(fp);
