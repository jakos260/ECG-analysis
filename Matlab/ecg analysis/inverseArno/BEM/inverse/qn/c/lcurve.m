echo off;
clear
hold off

fin=fopen('lcurve.tab','r');
nrow=fscanf(fin,'%d',[1,1]);
ncol=fscanf(fin,'%d',[1,1]);
datamat=fscanf(fin,'%f',[ncol,nrow]);
fclose(fin);
datamat=datamat';

loglog(datamat(:,3),datamat(:,7),'y-',datamat(:,3),datamat(:,7),'yo');
xlabel('reg');
ylabel('res');
%loglog(datamat(:,2),datamat(:,3),'y-',datamat(:,2),datamat(:,3),'yo');
%ylabel('reg');
%xlabel('res');
hold on
%axis([2.0 50.0 0.04 0.4 ]); 
%set(gca,'ytick',[0.05 0.1 0.2 0.3])
%set(gca,'xtick',[2 3 5 10 20 30 50])
%title('unrot L-curve udl');
%hold off

%orient landscape;
orient tall;
print lcurve.ps -deps





