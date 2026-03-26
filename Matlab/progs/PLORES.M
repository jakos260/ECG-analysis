echo off
fin=fopen('sinvals.ama','r');
nrow=fscanf(fin,'%d',[1,1]);
ncol=fscanf(fin,'%d',[1,1]);
a=fscanf(fin,'%f',[ncol,nrow]);
fclose(fin);
a=a';
t=(1:ncol);
plot(t(20:ncol),a(1,20:ncol))