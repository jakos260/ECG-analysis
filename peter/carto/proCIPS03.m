XML=xml2struct('C:\Users\peter\CIPS\ProCIPS03.xml');
wct=[557 290 302];
[Vv,I]=loadtri('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\mriModel\mriModel_ventricles.tri');
type=loadmat('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\mriModel\mriModel_ventricles.typ');

AMA=loadmat('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\mriModel\mriModel_thorax.vedl'); 
AMAV=40*AMA(end-length(Vv)+1:end,:);
AMAV= AMAV - ones(size(AMAV,1),1) * mean(AMA(wct,:));

dep = str2num(XML.CIPS{2}.CIPSDATA.ECG_cases.CIPS_ECG_CASE.ECGfiles.ECGFiledata{1}.selectedVentricularBeats.CIPS_RESULT_CASE.finaldep.Text);
t=0:round(max(dep))+10;
T= ones(length(dep),1)*t;   
SPECS.depSlope=1;
S=getSmode(T,dep,dep+300,SPECS,1);
EGM = AMAV *S;
focus = find(dep==min(dep));
qtriplot('delete *')

amplitude = max(EGM')-min(EGM');
qtriplot(Vv,I);
qtriplot(amplitude')


qtriplot('delete *')

load('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\Carto\Export\data\1-PVC RV_car.mat')
DATA.CAR.VER=DATA.CAR.VER(:,[3 1 2]);

VV= DATA.CAR.VER ;%.* str2num(D{1}.Camera.Aspectratio);
dd= Vv(410,:)-VV(641,:);

qtriplot(Vv,I)

qtriplot('trans 0.5')
qtriplot('linewidth 2')
load('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\Carto\Export\data\10-AO FAM_car.mat')

DATA.CAR.VER=DATA.CAR.VER(:,[3 1 2]);
VV=DATA.CAR.VER;
VV= [(VV(:,1)+dd(1)) (VV(:,2)+dd(2)) (VV(:,3)+dd(3))];
qtriplot(VV,DATA.CAR.ITRI)


load('C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS03\Carto\Export\data\1-PVC RV_car.mat')
DATA.CAR.VER=DATA.CAR.VER(:,[3 1 2]);


VV= DATA.CAR.VER ;%.* str2num(D{1}.Camera.Aspectratio);
dd= Vv(410,:)-VV(641,:);
VV= [(VV(:,1)+dd(1)) (VV(:,2)+dd(2)) (VV(:,3)+dd(3))];

% DATA.CAR.ITRI=DATA.CAR.ITRI(:,[1 3 2]);
[V,I]=make_sphere(2,0);
V= V*3.0;
a=find(DATA.TagPoints.Id == 9);
for i=1:length(a)
    b= VV(DATA.CAR.mapping(a(i),1),:);
    qtriplot([(V(:,1)+b(1)) (V(:,2) +b(2)) (V(:,3)+b(3))],I);
    qtriplot(['color ' num2str([1 0.2 0.2])])
end


qtriplot(VV,DATA.CAR.ITRI)
qtriplot(DATA.CAR.T * DATA.CAR.data(:,4) )
% qtriplot(VV,DATA.CAR.ITRI)
% VV= VV .* str2num(D{1}.Camera.Scale);
% qtriplot(DATA.CAR.VER,DATA.CAR.ITRI)


[A,D]=graphdist(DATA.CAR.ITRI,DATA.CAR.VER,4);
closetoABL = 206;
closetoABL = find(D(closetoABL,:) == min(D(closetoABL,DATA.CAR.mapping(:,1))));
% qtriplot(D(:,closetoABL))

dd= [(1:length(DATA.CAR.mapping))' D(DATA.CAR.mapping(:,1),closetoABL)];
dd=sortrows(dd,2);
% dd(dd(:,2)>0 & dd(:,2)<10,:)=[];


a=dd([1:4 6],1);
savetri('tmpel.tri',DATA.CAR.VER(DATA.CAR.mapping(a(2:end),1),:),[]);
qtriplot('file tmpel.tri')
figure(10)
clf
plot(EGM(focus,:),'r','linewidth',2);
hold on
i1=DATA.CAR.mapping(a(1),2);
plot(baselinecor(DATA.CAR.ECGM1(i1,1891:end-400)),'b','linewidth',2);
box off;grid on
legend('simulated (CIPS)', 'measured (Carto)')
ylabel('mV')
xlabel('ms')
title('EGM at Carto and CIPS estimated PVC location')


cols=['b  ';'k  ';'r  ';'g  ';'m  ';'y  ';'r: ';'k--';'g--';'c--';'m--';'y--'];
figure(1)
clf
leg=[];
for i=1:length(a)
    subplot(4,1,1)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('tip [mV]')
%     plot(DATA.CAR.ECGuni(i1,:),cols(i,:));
    hold on

    subplot(4,1,2)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('ring [mV]')
%     plot(DATA.CAR.ECGuni(i1,:),cols(i,:));
    hold on
    
    subplot(4,1,3)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)- DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('bipolar [mV]')
%     plot(DATA.CAR.ECGuni(i1,:),cols(i,:));
    hold on
    
    
    subplot(4,1,4)
    ECG=cell2mat(DATA.CAR.ECG_all(i1));
    plot(baselinecor(ECG(5,1850:end-400)),cols(i,:),'linewidth',2);
%     plot(ECG(5,:),cols(i,:));
    ylabel('V1 [mV]')
    hold on

    leg{i} = [num2str( round(dd(i,2)*10)/10 ) ' mm'];
end
subplot(4,1,1);grid
subplot(4,1,2);grid
subplot(4,1,3);grid
subplot(4,1,4);grid
legend(leg,'Location','southeast')


saveas(figure(1),'.\results\proCIPS_atInitmm.png')

dd(dd(:,2)>0 & dd(:,2)<10,:)=[];
% a=dd(1:6,1);
a=dd([1:4 ],1);

savetri('tmpel1.tri',DATA.CAR.VER(DATA.CAR.mapping(a,1),:),[]);
qtriplot('file tmpel1.tri')
qtriplot('color black')

figure(2)
leg=[];
for i=1:length(a)
    subplot(4,1,1)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('tip [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on

    subplot(4,1,2)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('ring [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on
    
    subplot(4,1,3)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)- DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('bipolar [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on
    
    
    subplot(4,1,4)
    ECG=cell2mat(DATA.CAR.ECG_all(i1));
    plot(baselinecor(ECG(5,1850:end-400)),cols(i,:),'linewidth',2);
%     plot(ECG(5,:),cols(i,:));
    ylabel('V1 [mV]')
    hold on

    leg{i} = [num2str( round(dd(i,2)*10)/10 ) ' mm'];
end
subplot(4,1,1);grid
subplot(4,1,2);grid
subplot(4,1,3);grid
subplot(4,1,4);grid
legend(leg,'Location','southeast')
saveas(figure(2),'.\results\proCIPS_at10mm.png')



dd(dd(:,2)>0 & dd(:,2)<20,:)=[];
a=dd([1 2 5],1);
savetri('tmpel2.tri',DATA.CAR.VER(DATA.CAR.mapping(a,1),:),[]);
qtriplot('file tmpel2.tri')
qtriplot('color red')
figure(3)
clf
leg=[];
for i=1:length(a)
    subplot(4,1,1)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('tip [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on

    subplot(4,1,2)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('ring [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on
    
    subplot(4,1,3)
    i1=DATA.CAR.mapping(a(i),2);
    plot(baselinecor(DATA.CAR.ECGM1(i1,1850:end-400)- DATA.CAR.ECGM2(i1,1850:end-400)),cols(i,:),'linewidth',2);
    ylabel('bipolar [mV]')
%     plot(DATA.CAR.ECGM1(i1,:),cols(i,:));
    hold on
    
    
    subplot(4,1,4)
    ECG=cell2mat(DATA.CAR.ECG_all(i1));
    plot(baselinecor(ECG(5,1850:end-400)),cols(i,:),'linewidth',2);
%     plot(ECG(5,:),cols(i,:));
    ylabel('V1 [mV]')
    hold on

    leg{i} = [num2str( round(dd(i,2)*10)/10 ) ' mm'];
end
subplot(4,1,1);grid
subplot(4,1,2);grid
subplot(4,1,3);grid
subplot(4,1,4);grid
legend(leg,'Location','southeast')
saveas(figure(3),'.\results\proCIPS_at20mm.png')

savetri('tmpel0.tri',DATA.CAR.VER(closetoABL,:),[]);
qtriplot('file tmpel0.tri')
qtriplot('color yellow')
