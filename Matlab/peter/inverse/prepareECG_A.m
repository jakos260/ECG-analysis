% Peter van Dam; 2013 July.
% select all ficucial points and estimate the dominatnt twave from these
% signals
% 20130915 oostep1: add Time_Vstim to specs, instruction in window title

function SPECS = prepareECG(BSM,LAY,varargin)

ECGextra = [];
doAtria = 0;
filename = [];
doSave = 1;
doCumsum =1;
doVstim=0;

if nargin < 2
	error('This routine needs at least two parameters');
else
	pp=1;
	while pp<=nargin-2
        key=lower(varargin{pp});
        switch key
            case 'atria'
                doAtria=varargin{pp+1};pp=pp+2;
            case 'extra'
                ECGextra=varargin{pp+1};pp=pp+2;
            case 'filename'
                filename=varargin{pp+1};pp=pp+2;
            case 'dosave'
                doSave=varargin{pp+1};pp=pp+2;
            case 'documsum'
                doCumsum = varargin{pp+1}; pp=pp+2;
            case 'dovstim'
                doVstim = varargin{pp+1}; pp=pp+2;
            otherwise
                error('unknown parameter');
        end
	end
end

% SUBJECT SPECIFIC INPUT SPECS
funtype=6; % product two logistic functions; one with shift
%funtype=7; % as type 6, but with added Gauss for U wave

if ~isempty(filename) && ~isempty(strfind(filename,'.spe') )
    fileout = [filename(1:strfind(filename,'.spe')-1) '.ecgspecs'];    
elseif ~isempty(filename)
    fileout = [filename '.ecgspecs'];
else
    fileout = [];
end
if exist(fileout,'file')
    saveFile = doSave;
else
    saveFile = 2;
end

if saveFile == 1 
    tmpSpecs = loadmat(fileout);
    SPECS.onsetP        = tmpSpecs(1);
    SPECS.onsetqrs      = tmpSpecs(2);
    SPECS.endtwave      = tmpSpecs(3);
    SPECS.time_Jpoint   = tmpSpecs(4);
    SPECS.time_apexT    = tmpSpecs(5);
    SPECS.time_apexU    = tmpSpecs(6);
    SPECS.endP          = tmpSpecs(7);
%     SPECS.initialSlope  = tmpSpecs(8);
%     SPECS.plateauslope  = tmpSpecs(9);
%     SPECS.repslope      = tmpSpecs(10);
%     SPECS.repCorrection = tmpSpecs(11);   
    SPECS.useCumsum     = tmpSpecs(12);   
    SPECS.qrsduration = SPECS.time_Jpoint;
    SPECS.qrstduration= SPECS.endtwave - SPECS.onsetqrs;
    if length(tmpSpecs) >=13
        SPECS.time_Vstim=tmpSpecs(13);
    else
        SPECS.time_Vstim=[];
    end
else
    SPECS = [];
end


% select teh fiducial point sin the ECG
SPECS = selectFiducucialPoints(SPECS, BSM,LAY, ECGextra,doAtria,doVstim);

SPECS.useCumsum = doCumsum;
SPECS.scaleAmpl     =[];
SPECS.qrsduration   = SPECS.time_Jpoint;
SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;

% compute source parameters
SPECS = estimatedFromTdom(SPECS,funtype);



if ~isempty(fileout)
    specs = zeros(11,1);
    specs(1) = SPECS.onsetP;
    specs(2) = SPECS.onsetqrs;
    specs(3) = SPECS.endtwave;
    specs(4) = SPECS.time_Jpoint;
    specs(5) = SPECS.time_apexT;
    specs(6) = SPECS.time_apexU;
    specs(7) = SPECS.endP;
%     specs(8) = SPECS.initialSlope ;
%     specs(9) = SPECS.plateauslope;
%     specs(10)= SPECS.repslope;
%     specs(11)= SPECS.repCorrection;   
    specs(12)= SPECS.useCumsum;
%     if length(specs) > 12
%         specs(13)=SPECS.time_Vstim;
%     end
    extra = 'onsetP onsetQRS endTwave Jpoint ApexT apexU depslope initialSlope plateauSlope repslope repCorrection useCumsum Vstim';
    saveasci(fileout,specs,extra)
end



%% ------------------------------------------------------------------------

function SPECS = selectFiducucialPoints(SPECS,ECGIN,LAYIN, ECGextra,doAtria,doVstim)
% inspect data
ECG = zeromean(ECGIN);
% ECG=baselinecor(subtracthum(ECGIN)); % oostep1
ECG(ECG>5)=2;
ECG(ECG<-5)=-2;


selectFiducial = isempty(SPECS);
fh=figure(301);
clf
if size(LAYIN,1) == 10
    LAY = [[1 9 0];[ones(9,1) (1:9)' (1:9)']];
    LAYUSE = LAY;
else
    LAY = LAYIN;
    nL = max(LAY(:,3));
    if nL < 30
        LAY = [[1 nL 0];[ones(nL,1) (1:nL)' (1:nL)']];
        LAYUSE =LAY;
    else
        nL = ceil(nL /3);
        LAY = [[1 nL 0];[ones(nL,1) (1:nL)' (1:3:nL*3)']];
        LAYUSE = LAY;%[[1 nL 0];[ones(nL,1) (1:nL)' (1:nL)']];
    end
end
if ~isempty(ECGextra)
    nL = LAY(1,2);
    nLext = size(ECGextra,1);
    LAYUSE(1,2) = nL + nLext;
    LAYUSE = [LAYUSE ; [ones(nLext,1) (nL+1:nL+nLext)' (nL+1:nL+nLext)']];
    ECG=[ECG;ECGextra];
end

sigscal = 10/max(max(abs(ECG)));
% sigplot(ECG,'',LAYUSE,sigscal,'b',1,0);
sigplot(ECG,'',LAYUSE,sigscal,'b');% oostep1 fix error on LAB. sigscal is amplification factor

hold on
% plot rmscurve for identifying/checking timing parameters
rrms = rms(ECG);
rrms = 2*rrms/min(1,max(rrms))*LAY(1,2);
nt = length(rrms);

t=(1:nt)/nt;
% plot(t,lowpassma(rrms,10),'r')
plot(t,rrms,'r'); % oostep1: not nice with pacing spikes

if selectFiducial && doAtria
    % identify major timing parameters QRST
    set(fh,'Name','onset Pwave?')
    disp('onset Pwave?')
    [x,~]=ginput(1);
    SPECS.onsetP=max(1,round(x/t(nt)*nt));
    
    set(fh,'Name','end Pwave?')
    disp('end Pwave?')
    [x,~]=ginput(1);
    SPECS.endP=max(SPECS.onsetP,round(x/t(nt)*nt));
    
elseif ~doAtria
    SPECS.onsetP = 1;
end
plot(t(SPECS.onsetP),rrms(SPECS.onsetP),'k*')
drawnow

if selectFiducial
    % identify major timing parameters QRST
    set(fh,'Name','onset QRS?')
    disp('onset QRS?')
    [x,~]=ginput(1);
    SPECS.onsetqrs=max(1,round(x/t(nt)*nt));
end
plot(t(SPECS.onsetqrs),rrms(SPECS.onsetqrs),'k*')


if selectFiducial
    set(fh,'Name','end T wave?')
    disp('end T wave?')
    [x,~]=ginput(1);
    SPECS.endtwave =min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.endtwave),rrms(SPECS.endtwave),'k*')


if doVstim && (selectFiducial || isempty(SPECS.time_Vstim))
    set(fh,'Name','Start Ventriclar stimulus? Right click to ignore');
    disp('Start Ventriclar stimulus? Right click to ignore');
    trms = rms(ECG);

    if max(trms) > 10*mean(trms)
        a = find(trms > 10*mean(trms));
        a(diff(a)<2)=[];
        if length(a) > 1
            SPECS.time_Vstim=a(2);
        else
            SPECS.time_Vstim=a(1);
        end
    else    
        [x,~,button]=ginput(1);
        if button==1
            SPECS.time_Vstim=max(1,round(x/t(nt)*nt));
        elseif button==3
            SPECS.time_Vstim=NaN;
        end
    end
end
if doVstim && ~isnan(SPECS.time_Vstim)
    plot(t(SPECS.time_Vstim),rrms(SPECS.time_Vstim),'r*');
end

ddrrms=diff(lowpassma(rrms,10),2);
sc=.5*max(rrms(SPECS.onsetqrs:SPECS.endtwave))/max(abs(ddrrms(max(1,SPECS.onsetqrs-2):SPECS.endtwave-2)));
plot(t(3:end),sc*ddrrms,'g');



PSI=baselinecor(ECG(:,SPECS.onsetqrs:SPECS.endtwave));
clf;
% sigplot(PSI,'',LAYUSE,sigscal,'b',1,0);
sigplot(PSI,'',LAYUSE,sigscal,'b');% oostep1 prevent LAB error in sigplot
hold on
rrms = rms(PSI);
showrrms = rrms * 2.0 /min(1,max(rrms))*LAY(1,2);

nt = size(PSI,2);
t= ( 1 : nt ) / nt;
plot(t,lowpassma(showrrms,5),'r')
% drrms = rms(diffrows(lowpassma(PSI,10)));
% plot(t,drrms* 2.0 /min(1,max(drrms))*LAY(1,2),'m')
hold on
if selectFiducial
    set(fh,'Name','end QRS?')
    disp('end QRS?')
    [x,~]=ginput(1);
    SPECS.time_Jpoint = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.time_Jpoint),showrrms(SPECS.time_Jpoint),'k*')
if selectFiducial
    set(fh,'Name','apex Twave?')
    disp('apex Twave?')
    [x,~]=ginput(1);
    SPECS.time_apexT = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.time_apexT), showrrms(SPECS.time_apexT),'k*')

if selectFiducial
    set(fh,'Name','apex Uwave?')
    disp('apex Uwave?')
    [x,~]=ginput(1);
    SPECS.time_apexU = min(round(x/t(nt)*nt),nt);
end
if SPECS.time_apexU > 0
    plot(t(SPECS.time_apexU), showrrms(SPECS.time_apexU),'k*')
end

t=(1:nt)/nt;
tdom = zeros(1,nt);
tdom(SPECS.time_Jpoint:end) = rrms(SPECS.time_Jpoint:end);
showtdom = 2*tdom/min(1,max(rrms))*LAY(1,2);
plot(t,showtdom ,'g')
plot([1 nt],[0 0],':k')

SPECS.tdom = tdom;
SPECS.rrms = rrms;


    


%% ------------------------------------------------------------------------
function SPECS = estimatedFromTdom(SPECS,funtype)

depSlope =2;

% initial estimates
pinit(1)= depSlope;     
pinit(2)= 0.01;
pinit(3)=-0.03;
pinit(4)= 0.06;
pinit(5)= SPECS.time_apexT - SPECS.time_Jpoint;   % all timing with respect to inset qrs
if funtype==6 
    pinit=pinit(1:5);
else
pinit(6)= SPECS.rrms(SPECS.time_apexU - SPECS.time_Jpoint);
pinit(7)= 50;  %std(uwave)
pinit(8)= SPECS.time_apexU - SPECS.time_Jpoint;
end

%tdom = baselinecor( lowpassma(SPECS.tdom(SPECS.time_Jpoint:SPECS.endtwave-SPECS.onsetqrs),20))';
tdom = lowpassma(SPECS.tdom(SPECS.time_Jpoint:SPECS.endtwave-SPECS.onsetqrs),20);
tdom = baselinecor(tdom,find(tdom== min(tdom(1:min(length(tdom),200)))),length(tdom))';


x    = (1:length(tdom))';
%% compute parameters for TMP from fit to dominant T wave
% signal from which tdom is constructed:
% fit analytical function (product of logistics) to tdomS
paramsCumsum = quamin_pvd(x,tdom,pinit,funtype);
paramsCumsum(1) = max(depSlope,paramsCumsum(1));
% paramsCumsum = quamin_pvd(x,tdom,paramsCumsum(2:5),funtype);
tdomFit = rfunc(x, paramsCumsum, 0, funtype);
tdomFit = max(SPECS.rrms(SPECS.time_Jpoint:end))/max(tdomFit)*tdomFit;
%tdomFit = baselinecor(tdomFit')';
figure(302)
clf
plot(x,tdom,'b')
hold on
plot(x,tdomFit,'k')


%% compute parameters for TMP from fit to dirved TMP waveform
APtwave = 1 - cumsum(tdom);
APtwave = APtwave - APtwave(end);
APtwave = APtwave / max(APtwave);

APmean = 1 - cumsum(tdomFit);
APmean = APmean - APmean(end);
APmean = APmean / max(APmean);
pinit(1) =  max(depSlope,paramsCumsum(1));
pinit(2) =  paramsCumsum(2);
pinit(3) =  paramsCumsum(3);
pinit(4) =  paramsCumsum(4);




if ~SPECS.useCumsum
    % fit analytical function (logistics) to tdom
    parms=quamin_pvd(x,APtwave,pinit,funtype);
    APfit = rfunc(x,parms,0,funtype);
    APfit = APfit - APfit(end);
    APfit = APfit / max(APfit);
    parms(5)     = parms(5)     - (SPECS.time_apexT - SPECS.time_Jpoint);
end
% paramsCumsum = [depSlope paramsCumsum];
paramsCumsum(5) = paramsCumsum(5) - (SPECS.time_apexT - SPECS.time_Jpoint);

if SPECS.useCumsum
    p = paramsCumsum;
else
    p = parms;
end

SPECS.depSlope      = p(1);
SPECS.initialSlope  = p(2);
SPECS.plateauslope  = p(3);
SPECS.repslope      = p(4);
SPECS.repCorrection = p(5);


dep = 10; 
rep = SPECS.time_apexT + dep;
t = 1:length(tdom)+100;



S  = getSmode(t, dep,    rep,       SPECS,4);
S1 = getSmode(t, dep,    rep - 100, SPECS,4);
S2 = getSmode(t, dep+100,rep - 100, SPECS,4);


% Savo  = getSmode(t, dep,rep,paramsCumsum,[],4,1);
% S1avo = getSmode(t, dep,rep-100,paramsCumsum,[],4,1);
% S2avo = getSmode(t, dep+100,rep-100,paramsCumsum,[],4,1);



figure(303);
clf;
subplot(2,1,1)
plot(APtwave,'b','linewidth',1)
hold on;
plot(APmean,'r','linewidth',1);

%plot(APfit,'k','linewidth',1);

axis tight
legend('1-int(Tdom)','1-int(fitted Tdom)')%,'fitted TMP')
subplot(2,1,2)
plot(t,S(1,:),'r','linewidth',1)
hold on
% plot(t,Savo(1,:),'k','linewidth',1)
% legend('direct getSmode','cumsum gets')
plot(t,S1(1,:),'r','linewidth',1)
plot(t,S2(1,:),'r','linewidth',1)
% plot(t,S1avo(1,:),'k','linewidth',1)
% plot(t,S2avo(1,:),'k','linewidth',1)
axis tight
%disp(num2str([paramsCumsum;parms]))


