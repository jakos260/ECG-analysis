% Peter van Dam; 2010 november.
% All rights reserved Peacs, Arnhem  the Netherlands
function prepare_geom_Haga(BSM,LAY,fileout)

% file prepare_geom.m
% date:071016

% view  and treat/prepare
% block 'specs' if new estimate is desired
% unblock 'fileout' if resulting data need to be stored

% clear all

% SUBJECT SPECIFIC INPUT SPECS

funtype=15;

PHI = BSM(1:size(LAY)-1,:);
PSI = PHI;

% inspect data
figure(301)
clf

sigscal=max(max(abs(PSI)))/2;
sigplot(PSI,'',LAY,sigscal,'b',1,0);
hold on
% plot rmscurve for identifying/checking timing parameters
rrms=rms(PSI);
rrms=2*rrms/max(rrms)*LAY(1,2);
nt=length(rrms);

t=(1:nt)/nt*LAY(1,1);
plot(t,rrms,'r')

% identify major timing parameters QRST

disp('onset QRS?')
[x,y]=ginput(1);
onsetqrs=max(1,round(x/t(nt)*nt));
plot(t(onsetqrs),rrms(onsetqrs),'k*')

disp('end T wave?')
[x,y]=ginput(1);
endtwave =min(round(x/t(nt)*nt),nt);
plot(t(endtwave),rrms(endtwave),'k*')


figure(302); clf
PSI=baselinecor(PSI,onsetqrs,endtwave);
rrms=rms(PSI);
rrms=rrms/max(rrms(onsetqrs+150:nt));
plot(rrms(onsetqrs:endtwave),'r')
hold on



disp('end QRS?')
[x,y]=ginput(1);
endqrs=round(x)+onsetqrs;
clf
plot(rrms,'r');hold on
plot(onsetqrs,rrms(onsetqrs),'k*')
plot(endqrs,rrms(endqrs),'k*')
plot(endtwave,rrms(endtwave),'k*')

% ttm: meanrep estimated from timing apex RMS
[ma,ttm]=max(rrms(endqrs:nt));
apext=ttm+endqrs-1;
plot(apext,rrms(apext),'*')
qrsdur=endqrs-onsetqrs+1;

%% compute parameters for TMP from fit to dominant T wave
% signal from which tdom is constructed:
ntt=endtwave-onsetqrs+1;
tdom=ones(1,ntt)*PSI(:,onsetqrs:endtwave)'*PSI(:,onsetqrs:endtwave);
tdom=tdom*rrms(apext)/max(tdom(qrsdur:length(tdom)));
% tdom=tdom/sum(tdom);
% plot(onsetqrs:endtwave,tdom,'m')
plot([1 nt],[0 0],':k')

tbeg=onsetqrs;% endqrs;
[tdom,ttm]=gettdompvd(PSI,onsetqrs, endqrs);
plot(1:length(tdom),tdom*max(rrms(endqrs:nt))/max(abs(tdom)),'g')
%
% tbeg=endqrs;%+ttm-(length(tdom)-ttm);
[tdom,ttm]=gettdompvd(PSI(:,tbeg:endtwave),1,qrsdur);
%
%
pinit(1)= 1;%tdom(1);	% overall scaling factor
pinit(2)= 1;%0.01;	% the initial value of phiup, phidown,phinotch
pinit(3)=.01;
pinit(4)=0.22;
pinit(5)=ttm;

APmean=[1-cumsum(tdom) zeros(1,100)]';
% fit analytical function (logistics) to tdom
parms=quamin_pvd((1:length(APmean))',APmean,pinit,15);
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
legend('1-int(Tdom)','fitted','getSmode')

%%
% The extra factors (1.2 and 2) were found after comparing the Vm of beat 7
% (non uniform cell types without Ito) beta 7 was provided by Mark Potse.
onsetP = size(PHI,2);
% specs=[sigscal onsetqrs endqrs ttm endtwave [parms(3) parms(4)] parms(5)-ttm ]';
specs=[sigscal onsetqrs endqrs ttm onsetP [parms(3) parms(4)] parms(5)-ttm ]';

% disp(num2str(specs',4))

saveasci(fileout,specs)




