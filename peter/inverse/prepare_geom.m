% Peter van Dam; 2013 July. 
function GEOM=prepare_geom(GEOM,fileout,saveFile,varargin)

% file prepare_geom.m
% date:071016

% view  and treat/prepare  
% block 'specs' if new estimate is desired
% unblock 'fileout' if resulting data need to be stored

% clear all

% SUBJECT SPECIFIC INPUT SPECS
funtype=6; % product two logistic functions; one with shift
%funtype=7; % as type 6, but with added Gauss for U wave

%default values
GEOM.specs(6)=0.045517;
GEOM.specs(7)=0.044909; 
GEOM.specs(8)=0.077413;
ECGextra = [];
if nargin > 3
    ECGextra = varargin{4};
end
if saveFile == 1 && exist(fileout,'file')
	GEOM.specs=loadmat(fileout);
	GEOM.pS=GEOM.specs(6:end);
    GEOM.specs(5) = min(GEOM.specs(5),size(GEOM.BSM,2));
    GEOM.qrsduration = GEOM.specs(3) - GEOM.specs(2)+1; % time_Jpoint - onsetqrs  
    GEOM.qrstduration=size(GEOM.BSM,2)-GEOM.specs(2)+1;
end
SPECS = selectFiducucialPoints(GEOM, ECGextra,saveFile);

PSI = GEOM.BSM(:,SPECS.onsetqrs:SPECS.endtwave);

% initial estimates
p(1)= 4;     % overall scaling factor T dom
p(2)= 0.01;
p(3)=-0.03;
p(4)= 0.04;
p(5)= SPECS.time_apexT;   % all timing with respect to inset qrs

p(6)= SPECS.rrms(SPECS.time_apexU);
p(7)= 50;  %std(uwave)
p(8)= SPECS.time_apexU;

y = SPECS.tdom(SPECS.time_Jpoint : length(SPECS.rrms))';
x =(SPECS.time_Jpoint : length(SPECS.rrms))';

pinit=p(1:5);
if funtype==7, 
    pinit=p(1:8); 
end
%% compute parameters for TMP from fit to dominant T wave   
% signal from which tdom is constructed:
% fit analytical function (product of logistics) to tdom
parms = quamin_pvd((1:length(SPECS.rrms))',SPECS.tdom',pinit,funtype);




x=1:length(SPECS.rrms);
y1=rfunc(x,parms,0,funtype);
y1=max(SPECS.rrms(SPECS.time_Jpoint:end))/max(y1)*y1;
plot(x,y1,'k')
tmp=1-cumsum(y1);
tmp=tmp-tmp(end);
tmp=tmp/max(tmp);



%% compute parameters for TMP from fit to dominant T wave   
APmean=[1-cumsum(SPECS.tdom) zeros(1,100)]';
% fit analytical function (logistics) to tdom
parms=quamin_pvd((1:length(APmean))',APmean,pinit,funtype);
% parms=quamin_pvd((1:length(tdom))',tdom',pinit,6)
yap=rfunc_pvd(1:length(APmean),parms,0,funtype);
figure(303);  set(gcf,'position',[63   363   485   400]); clf;
plot(APmean,'b');hold on;plot(yap,'k');
dep=10; rep=ttm+dep+(parms(5)-ttm);
tmp=[zeros(1,dep) 1-cumsum(tdom) zeros(1,nt-endtwave)]';

if funtype==12
	t=1:length(tmp);S=getSvdeprep(t,dep,rep,[parms(3) parms(4)],0,[],3);
else
t=1:length(tmp);S=getSmode(t,dep,rep,[parms(3) parms(4)],[],4);
end
plot(t,S(1,:),'r')
axis tight
title(['subject: ' GEOM.subject]);legend('1-int(Tdom)','fitted','getSmode')  

%%
% The extra factors (1.2 and 2) were found after comparing the Vm of beat 7
% (non uniform cell types without Ito) beta 7 was provided by Mark Potse.  
onsetP = size(PHI,2);
% specs=[sigscal onsetqrs time_Jpoint ttm endtwave [parms(3) parms(4)] parms(5)-ttm ]';
specs=[sigscal onsetqrs time_Jpoint ttm+tbeg onsetP [parms(3) parms(4)] parms(5)-ttm ]';
GEOM.specs = specs;    
GEOM.pS = specs(6:end);
GEOM.qrsduration=GEOM.specs(3)-GEOM.specs(2)+1; % time_Jpoint - onsetqrs  
GEOM.qrstduration=size(GEOM.BSM,2)-GEOM.specs(2)+1;

disp(num2str(specs',4))
if saveFile
%     resp=questdlg('store downslope parameters?','dominat T','yes','no','no');
%     if strcmp(resp,'yes')
        saveasci(fileout,specs)

%     end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SPECS = selectFiducucialPoints(GEOM, ECGextra,saveFile)
% inspect data
ECG = GEOM.BSM;
LAYIN = GEOM.LAY;

figure(301)
clf
if size(LAYIN,1) == 10
    LAY = [[1 9 0];[ones(9,1) [1:9]' [1:9]']];
    LAYUSE = LAY;
else
    LAY = LAYIN;
    nL = max(LAY(:,3));
    if nL < 30
        LAY = [[1 nL 0];[ones(nL,1) [1:nL]' [1:nL]']];
        LAYUSE =LAY;
    else
        nL = ceil(nL /3);
        LAY = [[1 nL 0];[ones(nL,1) [1:nL]' [1:3:nL*3]']];
        LAYUSE = [[1 nL 0];[ones(nL,1) [1:nL]' [1:nL]']];
    end
end
if ~isempty(ECGextra)
    nL = LAY(1,2);
    nLext = size(ECGextra,1);
    LAYUSE(1,2) = nL + nLext;
    LAYUSE = [LAYUSE ; [ones(nLext,1) (nL+1:nL+nLext)' (nL+1:nL+nLext)']];
    ECG=[ECG;ECGextra];
end

sigscal=max(max(abs(ECG))) * 2;
sigplot(ECG,'',LAYUSE,sigscal,'b',1,0);
hold on
% plot rmscurve for identifying/checking timing parameters
rrms = rms(ECG);
rrms = 2*rrms/min(1,max(rrms))*LAY(1,2);
nt = length(rrms);

t=(1:nt)/nt;
plot(t,rrms,'r')
if saveFile
    % identify major timing parameters QRST
    disp('onset QRS?')
    [x,~]=ginput(1);
    SPECS.onsetqrs=max(1,round(x/t(nt)*nt));
end
    
plot(t(SPECS.onsetqrs),rrms(SPECS.onsetqrs),'k*')
if saveFile
    disp('end T wave?')
    [x,~]=ginput(1);
    SPECS.endtwave =min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.endtwave),rrms(SPECS.endtwave),'k*')

% figure(302); clf
clf
sigplot(ECG,'',LAYUSE,sigscal,'b',1,0);
hold on

PSI=baselinecor(ECG(:,SPECS.onsetqrs:SPECS.endtwave));
rrms = rms(PSI);
nt = size(PSI,2);
t=(1:nt)/nt;
plot(t,2*rrms/min(1,max(rrms))*LAY(1,2),'r')
hold on       
if saveFile
    disp('end QRS?')
    [x,~]=ginput(1);
    SPECS.time_Jpoint = min(round(x/t(nt)*nt),nt);
end
plot(SPECS.time_Jpoint,rrms(SPECS.time_Jpoint),'k*')
if saveFile
    disp('apex Twave?')
    [x,~]=ginput(1);
    SPECS.time_apexT = min(round(x/t(nt)*nt),nt);
end
plot(SPECS.time_apexT,rrms(SPECS.time_apexT),'k*')

if saveFile
    disp('apex Uwave?')
    [x,~]=ginput(1);
    SPECS.time_apexU = min(round(x/t(nt)*nt),nt);
if saveFile
plot(SPECS.time_apexU,rrms(SPECS.time_apexU),'k*')

t=(1:nt)/nt;
tdom = zeros(1,nt);
tdom(SPECS.time_Jpoint:end) = rrms(SPECS.time_Jpoint:end);
plot(t,2*tdom/min(1,max(rrms))*LAY(1,2),'g')


plot([1 nt],[0 0],':k')

SPECS.tdom = tdom;
SPECS.rrms = rrms;



% clf
% plot(rrms,'r');hold on
% % plot(onsetqrs,rrms(onsetqrs),'k*')
% plot(SPECS.onsetqrs,rrms(1),'k*')
% plot(SPECS.time_Jpoint,rrms(SPECS.time_Jpoint-SPECS.onsetqrs),'k*')
% plot(SPECS.endtwave,rrms(SPECS.endtwave - SPECS.onsetqrs),'k*')
%     
% % ttm: meanrep estimated from timing apex RMS 
% % [~,ttm]=max(rrms(time_Jpoint:nt));
% [~,ttm]=max(rrms);
% SPECS.apext=ttm;%+time_Jpoint-1;
% SPECS.qrsduration = SPECS.time_Jpoint;


