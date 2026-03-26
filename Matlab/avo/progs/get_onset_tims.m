% get_onset_tims.m
% function TIMS=get_tims(PHI,icase,t,tbeg,tend,rmsmode,mode);
% use unweighted, common reference data  only
% 2013_01_06
function TIMS=get_tims(PHI,icase,t,tbeg,tend,rmsmode,mode);
global sint
sampling_rate=round(1000/sint)

nt=size(PHI,2);
t=1:nt;

% figure(2)
% clf
% start identifying timing markers onset QRS; unaugmented extremities assumed

PSI=PHI;  

single_1=rms(PSI,rmsmode); 

'STD-curve of unfiltered data'

high_pass=5;  

winh=sampling_rate/high_pass  % uniform moving average window
HPS=lowpassma(PSI,winh);
HPS=PSI-HPS;
['highpass frequency: ' num2str(high_pass) 'Hz']
single_2=rms(HPS,rmsmode);

low_pass=50;
winl=sampling_rate/low_pass % uniform moving average window
['lowpass frequency: ' num2str(low_pass) 'Hz']
single_2=lowpassma(single_2,winl);

maxsingle1=max(abs(single_1(tbeg:tend)));
maxsingle2=max(abs(single_2(tbeg:tend)));
maxsingle=max(maxsingle1,maxsingle2);
scal=1/maxsingle;

figure(2)
clf
singleplot_2=plot(tbeg:tend,scal*single_2(tbeg:tend),'b');
title(['STD_curve of bandpass-filtered signal; icase: ' num2str(icase)],'interpreter','none')
hold on

PEAKS=[];
% select a threshold, based on the first 2 seconds
threshold=0.5*max(single_2(1:round(2000/sint)));
PEAKS=   getpeaks(single_2(1:round(2000/sint)),threshold,80/sint);
% PEAKS selects first of extremes in intervals of 80 ms around highest local extremes

if isempty(PEAKS),
    beep;
    'peak detection failed; check signal quality and/or threshold settings'
    pause
else
    PEAKS=getpeaks(single_2,threshold,80/sint);
end

if ~isempty(PEAKS),
    % find crude estimate for onset QRS
    nk=size(PEAKS,1);
    npeaks=nk
    k=1:nk;
    onsets=[];
    for i=1:nk,
        est(i)       =PEAKS(k(i),1);
        onsets(i)    =PEAKS(k(i),2)-1;
        while single_2(onsets(i))<=est(i) | single_2(onsets(i))>PEAKS(k(i),1)/10,
            est(i)   = single_2(onsets(i));
            if onsets(i)==1, break,end
            onsets(i)= onsets(i)-1;
        end
        if onsets(i)>tbeg & onsets(i)<=tend,
            plot(onsets(i),scal*est(i),'k+','linewidth',1.5)
        end
    end
   
    zerotims=onsets';
    zerotims(zerotims<=1)=[];
    zerotims=unique(zerotims);
    nz=size(zerotims,1);
    
    intervals=zerotims(2:end)-zerotims(1:end-1);
    std_intervals=std(intervals,1)*sint;
    mean_intervals=mean(intervals)*sint;
    ['intervals mean and std: ' num2str(mean_intervals) ' ; ' num2str(std_intervals) '  (ms)']
    if std_intervals<20,
        beep
        'low std_intervals may indicate pacemaker'
        pause
    end
end

pause

TEST=PHI;
TEST=lowpassma(TEST,winl);
TEST=baselsplines(TEST,zerotims,2);
testsignal=rms(TEST,rmsmode);

nshow=tend-tbeg+1

itbeg=1;

while itbeg<nt,
    itend=min(itbeg+nshow-1,nt);
    figure(2)
    clf
    plot(itbeg:itend,scal*testsignal(itbeg:itend),'b');
    hold on
    izt =zerotims(zerotims<=itend&zerotims>=itbeg);
    plot(izt,scal*testsignal(izt),'k*');
    title(['STD_curve lowpass-filtered signal;crude-BLC; icase: ' num2str(icase)],'interpreter','none')
    
    deriv=diffrows(testsignal);
    threshold=0.5*max(deriv);
    PEAKS2=getpeaks(deriv,threshold,50/sint,-1);
    
    upslopetims=PEAKS2(:,2);
    upslopetims=unique(upslopetims);
    %plot(iupt,scal*testsignal(iupt),'m*');
   
    if upslopetims(1)<zerotims(1), upslopetims(1)=[]; end
    if ismember(zerotims(end),upslopetims(end-1):upslopetims(end))==0, upslopetims(end)=[];end

    iupt=upslopetims(upslopetims<=itend&upslopetims>itbeg);
    plot(iupt,scal*testsignal(iupt),'m*');
    izt =zerotims(zerotims<=itend&zerotims>itbeg);
    hold on
    plot(izt,scal*testsignal(izt),'k*');
    
    ipt=PEAKS(PEAKS(:,2)<=itend&PEAKS(:,2)>=itbeg,2);
    plot(ipt,scal*testsignal(ipt),'*r')
   
    tend
    'check that each zero marker (black) is followed by an upslope marker'
    itbeg=itend+1;
    nzeros=size(izt,1)
    nups=size(iupt,1)
    pause
end

nz=size(zerotims,1)
nup=size(upslopetims,1)

if nup~=nz,
    beep
 'nup should equal nz'   
pause
end

TIMS=[ones(nz,1) zerotims  upslopetims ];

for k=1:nz;
    DIPS=[];
    itest=0;
    low=inf;
    
    while itest<upslopetims(k)-zerotims(k),
        itest=itest+1;
        ztims=zerotims;
        ztims(k)=upslopetims(k)-itest;
        TEST=lowpassma(PHI,winl);
        
        TEST=baselsplines(TEST,ztims,2);
        testsignal=rms(TEST,rmsmode);
        
        dips=scal*testsignal(ztims-1)+testsignal(ztims+1);
        DIPS=[DIPS dips'];
        
        if dips(k)< low,
            low=dips(k);
        else
            break,
        end 
    end
    zerotims(k)=ztims(k);
end

nz=size(zerotims,1);
figure(2)
clf
plot(tbeg:tend,scal*testsignal(tbeg:tend))
hold on

izt =zerotims(zerotims<=itend&zerotims>itbeg);
plot(izt,scal*testsignal(izt),'k*');
 
title(['STD_curve of lowpass-, BLCorrected signals; icase: ' num2str(icase)],'interpreter','none')

intervals=zerotims(2:end)-zerotims(1:end-1);
std_intervals=std(intervals,1)*sint;
mean_intervals=mean(intervals)*sint;
['intervals mean and std: ' num2str(mean_intervals) ' ; ' num2str(std_intervals) '  (ms)']

pause

TIMS=zerotims;

