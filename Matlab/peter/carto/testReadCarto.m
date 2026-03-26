clear

basepath ='C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS'

pathname='proCIPS02\carto\Patient 2014_08_18\Study 1\Export';xmlin='Study 1 09_02_2014 07-57-09.xml';
pathname='proCIPS03\Carto\Export'; xmlin='Study 1 09_02_2014 07-47-07.xml';
pathname='proCIPS17\Carto\Export21April2015';xmlin='Study 1 05_05_2015 07-23-03.xml';
D = readCartoExportData(fullfile(basepath,pathname), xmlin);








return
clf
subplot(2,1,1)
ii1=576;
ii2=1345;
i1=find(DATA.mapping(:,1)==ii1);
i2=find(DATA.mapping(:,1)==ii2);

plot(DATA.ECGuni(i1,:));
hold on
plot(DATA.ECGuni(i2,:),'r');

subplot(2,1,2)
ECG=cell2mat(DATA.ECG_all(i1));
plot(ECG(5,:));
hold on
ECG=cell2mat(DATA.ECG_all(i2));
plot(ECG(5,:),'r');


savetri('tmpel.tri',DATA.VER([ii1 ii2],:),[]);

qtriplot('delete *');
qtriplot(DATA.VER,DATA.ITRI);
% qtriplot(['marker vertex ' num2str([ ii1 ii2])])
A= zeros(size(ECG));
A(DATA.mapping(:,1),:)=DATA.ECGuni(:,1880:2100);
ECG = DATA.T * DATA.ECGuni(:,1880:2100);
qtriplot(ECG)
qtriplot('file tmpel.tri')
for i=1:size(DATA.data,2)
    B= zeros(size(DATA.VER,1),1);
    B(DATA.mapping(:,1),:) = DATA.data(:,i);
    qtriplot(B)
    range(B)
%     qtriplot(DATA.T * DATA.data(:,i))
    qtriplot('funscale auto')
    pause
end