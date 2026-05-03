% function SPECS = prepareECG(BSM,LAY,varargin)
%
% Determine the parameters that describe the TMP function form.
% NOT the actual depolarization and repolarization times!!
%
% INPUT:
% BSM           = BSM
% LAY           = layoutfile
%
% OPTIONAL INPUT:
% 'atria'       = [1] P-wave included in analyses for atria [default = 0] 
% 'extra'       = additional ECG to be included in analyses [default = [] ]
% 'filename'	= filename of bsm-beat, without extension! [default = [] ]
% 'dosave'      = [2] rewrite ecgspecs, [1] read ecgspecs when available [default = 1]
% 'docusum'     = [0] fit cumulative rms Twave; better, but slower, then default [default = 1]
% 'dovstim'     = [1] Option to indicate start ventricular stimulus [default = 0]
% 'funtype'     = [6] product two logistic functions; one with shift. [7] same as type 6, but with added Gauss for U wave
%
% Peter van Dam; 2013 July.
% select all ficucial points and estimate the dominant twave from these signals
%
% 20130915 oostep1: add Time_Vstim to specs, instruction in window title
%
% AJ: Feb. 2, 2015: added peak QRS code.
%

function SPECS = prepareECG(BSM,LAY,varargin)

%% 
% default input:
ECGextra    = [];
doAtria     = 0;
filename    = [];
doSave      = 1;
doCumsum    = 0;
doVstim     = 0;
funtype     = 6; % [6] product two logistic functions; one with shift. [7] as type 6, but with added Gauss for U wave

% check input and replace default:
if nargin < 2,
    error('This routine needs at least two parameters');
else
    pp = 1;
    while pp <= nargin-2
        key = lower(varargin{pp});
        switch key
            case 'atria'
                doAtria = varargin{pp+1};   pp=pp+2;
            case 'extra'
                ECGextra = varargin{pp+1};  pp=pp+2;
            case 'filename'
                filename = varargin{pp+1};  pp=pp+2;
            case 'dosave'
                doSave = varargin{pp+1};    pp=pp+2;
            case 'documsum'
                doCumsum = varargin{pp+1};  pp=pp+2;
            case 'dovstim'
                doVstim = varargin{pp+1};   pp=pp+2;
            case 'funtype'
                funtype = varargin{pp+1};   pp=pp+2;
            otherwise
                error('unknown parameter');
        end
    end
end

% create output filename:
if ~isempty(filename) && ~isempty(strfind(filename,'.spe')),
    fileout = [filename(1:strfind(filename,'.spe')-1) '.ecgspecs'];
elseif ~isempty(filename)
    fileout = [filename '.ecgspecs'];
else
    fileout = [];
end

% check if *.ecgspecs exist, end set if it has to be overwritten or not!
if exist(fileout,'file'), saveFile = doSave; else saveFile = 2; end

% when *.ecgspecs exists and has to be used, load SPECS.
if saveFile == 1,
    tmpSpecs            = loadmat(fileout);
    SPECS.onsetP        = tmpSpecs(1);
    SPECS.onsetqrs      = tmpSpecs(2);
    SPECS.endtwave      = tmpSpecs(3);
    SPECS.time_Jpoint   = tmpSpecs(4);
    SPECS.time_apexT    = tmpSpecs(5);
    SPECS.time_apexU    = tmpSpecs(6);
    SPECS.useCumsum     = tmpSpecs(12);
    SPECS.qrsduration   = SPECS.time_Jpoint;
    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
    if length(tmpSpecs) >= 13, SPECS.time_Vstim = tmpSpecs(13); else SPECS.time_Vstim = []; end
    if length(tmpSpecs) >= 14, SPECS.peakqrs    = tmpSpecs(14); else SPECS.peakqrs    = []; end   % added peak QRS 20150202
else
    SPECS = [];
end

%%
% select the fiducial points in the ECG
SPECS               = selectFiducucialPoints(SPECS,BSM,LAY,ECGextra,doAtria,doVstim);   % indication of points in timesignal [rms BSM]
SPECS.useCumsum     = doCumsum;
SPECS.scaleAmpl     = [];
SPECS.qrsduration   = SPECS.time_Jpoint;                                                % duration QRS
SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;                                  % duration QRST

% compute source parameters
SPECS = estimatedFromTdom(SPECS,funtype);

% save results to fileout:
if ~isempty(fileout),
    specs       = zeros(11,1);
    specs(1)    = SPECS.onsetP;
    specs(2)    = SPECS.onsetqrs;
    specs(3)    = SPECS.endtwave;
    specs(4)    = SPECS.time_Jpoint;
    specs(5)    = SPECS.time_apexT;
    specs(6)    = SPECS.time_apexU;
    specs(7)    = SPECS.depSlope;
    specs(8)    = SPECS.initialSlope ;
    specs(9)    = SPECS.plateauslope;
    specs(10)   = SPECS.repslope;
    specs(11)   = SPECS.repCorrection;
    specs(12)   = SPECS.useCumsum;
    if ~isempty(SPECS.time_Vstim), specs(13) = SPECS.time_Vstim;    else specs(13) = nan; end
    if ~isempty(SPECS.peakqrs),    specs(14) = SPECS.peakqrs;       else specs(14) = nan; end % added peak QRS 20150202
    extra = 'onsetP onsetQRS endTwave Jpoint ApexT apexU depslope initialSlope plateauSlope repslope repCorrection useCumsum Vstim peakQRS'; % added peak QRS 20150202
    saveasci(fileout,specs,extra)
end
 
% ------------------------------------------------------------------------ %

function SPECS = selectFiducucialPoints(SPECS,ECGIN,LAYIN, ECGextra,doAtria,doVstim)

% selectFiducucialPoints(SPECS,BSM,LAY,ECGextra,doAtria,doVstim)

% inspect data
ECG             = zeromean(ECGIN);  % ECGIN = BSM
ECG(ECG>5)      = 2;                % ECG value above 5, make 2 [cutt-off]
ECG(ECG<-5)     = -2;               % ECG value under -5, make -2 [cutt-off]
selectFiducial  = isempty(SPECS);   % check if SPECS are loaded or not
fh              = figure(301);
clf

% make visualization for BSM, based on layout input
if size(LAYIN,1) == 10,
    LAY     = [[1 9 0];[ones(9,1) (1:9)' (1:9)']];
    LAYUSE  = LAY;
else
    LAY = LAYIN;
    nL  = max(LAY(:,3));
    if nL < 30,
        LAY     = [[1 nL 0];[ones(nL,1) (1:nL)' (1:nL)']];
        LAYUSE  = LAY;
    else
        nL      = floor(nL /3);
        LAY     = [[1 nL 0];[ones(nL,1) (1:nL)' (1:3:nL*3)']];
        LAYUSE  = LAY;
    end
end

% if ECGextra is not empty, add ECGextra to ECG and layout
if ~isempty(ECGextra)
    nL          = LAY(1,2);
    nLext       = size(ECGextra,1);
    LAYUSE(1,2) = nL + nLext;
    LAYUSE      = [LAYUSE ; [ones(nLext,1) (nL+1:nL+nLext)' (nL+1:nL+nLext)']];
    ECG         = [ECG;ECGextra];
end

% visualize ECG, using layout.
sigscal = 10/max(max(abs(ECG)));
sigplot(ECG,'',LAYUSE,sigscal,'b'); % oostep1 fix error on LAB. sigscal is amplification factor
hold on

% plot rmscurve for identifying/checking timing parameters
rrms    = rms(ECG);
rrms    = 2*rrms/min(1,max(rrms))*LAY(1,2);
nt      = length(rrms);

t = (1:nt)/nt;
plot(t,rrms,'r'); % oostep1: not nice with pacing spikes

if selectFiducial && doAtria,
    % identify major timing parameters QRST
    set(fh,'Name','onset Pwave?')
    disp('onset Pwave?')
    [x,ytemp]       = ginput(1);
    SPECS.onsetP    = max(1,round(x/t(nt)*nt));
else
    SPECS.onsetP = 1;
end
plot(t(SPECS.onsetP),rrms(SPECS.onsetP),'k*')
drawnow

if doVstim && (selectFiducial || isempty(SPECS.time_Vstim)),
    set(fh,'Name','Start Ventriclar stimulus? Right click to ignore. Any key accepts preset value.');
    disp('Start Ventriclar stimulus? Right click to ignore. Any key accepts preset value.');
    trms                = rms(ECG);
    SPECS.time_Vstim    = NaN;
    
    if max(trms) > 10*mean(trms),
        a                       = find(trms > 10*mean(trms));
        a(1+find(diff(a)<10))   = [];
        if length(a) > 1,
            SPECS.time_Vstim = a(2)-1;
        else
            SPECS.time_Vstim = a(1)-1;
        end
        plot(t(SPECS.time_Vstim),min(3,rrms(SPECS.time_Vstim)),'r*');
    end
    
    [x,ytemp,button] = ginput(1);
    
    % check button input [left or right mouseclick]
    if button == 1,
        SPECS.time_Vstim = max(1,round(x/t(nt)*nt));
    elseif button == 3,
        SPECS.time_Vstim = NaN;
    end
    
end

if selectFiducial,
    % identify major timing parameters QRST
    set(fh,'Name','onset QRS?')
    disp('onset QRS?')
    [x,ytemp] = ginput(1);
    SPECS.onsetqrs = max(1,round(x/t(nt)*nt));
end
plot(t(SPECS.onsetqrs),rrms(SPECS.onsetqrs),'k*')

if selectFiducial,
    set(fh,'Name','end T wave?')
    disp('end T wave?')
    [x,ytemp] = ginput(1);
    SPECS.endtwave = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.endtwave),rrms(SPECS.endtwave),'k*')

if doVstim && ~isnan(SPECS.time_Vstim),
    plot(t(SPECS.time_Vstim),rrms(SPECS.time_Vstim),'r*');
end

ddrrms = diff(lowpassma(rrms,10),2);
if SPECS.onsetqrs > 2, adj_onset = SPECS.onsetqrs-2; else adj_onset = 1; end
sc = .5*max(rrms(SPECS.onsetqrs:SPECS.endtwave))/max(abs(ddrrms(adj_onset:SPECS.endtwave-2)));
plot(t(3:end),sc*ddrrms,'g');

PSI = baselinecor(ECG(:,SPECS.onsetqrs:SPECS.endtwave));
clf;

sigplot(PSI,'',LAYUSE,sigscal,'b');                 % oostep1 prevent LAB error in sigplot
hold on
rrms        = rms(PSI);                             % rms of baseline corrected signal [ECG = BSM]
showrrms    = rrms*2.0 /min(1,max(rrms))*LAY(1,2);  % rms for visualization

nt  = size(PSI,2);
t   = ( 1 : nt ) / nt;
plot(t,lowpassma(showrrms,5),'r')

hold on
if selectFiducial,
    set(fh,'Name','end QRS?')
    disp('end QRS?')
    [x,ytemp]           = ginput(1);
    SPECS.time_Jpoint   = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.time_Jpoint),showrrms(SPECS.time_Jpoint),'k*')

if selectFiducial,
    set(fh,'Name','apex Twave?') % peak Twave
    disp('apex Twave?')
    [x,ytemp]           = ginput(1);
    SPECS.time_apexT    = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.time_apexT), showrrms(SPECS.time_apexT),'k*')

if selectFiducial,
    set(fh,'Name','apex Uwave?') % Uwave, after T
    disp('apex Uwave?')
    [x,ytemp]           = ginput(1);
    SPECS.time_apexU    = min(round(x/t(nt)*nt),nt);
end
plot(t(SPECS.time_apexU), showrrms(SPECS.time_apexU),'k*')

% added peak QRS 20150202
if selectFiducial,
    set(fh,'Name','peak QRS?') % Peak QRS
    disp('peak QRS?')
    [x,ytemp]           = ginput(1);
    SPECS.peakqrs       = min(round(x/t(nt)*nt),nt);
end
if ~isempty(SPECS.peakqrs) && ~isnan(SPECS.peakqrs), plot(t(SPECS.peakqrs), showrrms(SPECS.peakqrs),'k*'); end

t                           = (1:nt)/nt;
tdom                        = zeros(1,nt);
tdom(SPECS.time_Jpoint:end) = rrms(SPECS.time_Jpoint:end);
showtdom                    = 2*tdom/min(1,max(rrms))*LAY(1,2);
plot(t,showtdom ,'g')
plot([1 nt],[0 0],':k')

SPECS.tdom = tdom;  % rms of signal from Jpoint till end
SPECS.rrms = rrms;  % rms of complete signal

%% ------------------------------------------------------------------------ %
function SPECS = estimatedFromTdom(SPECS,funtype)

% initial estimates: pinit = [depslope, initial value Tdom, (pos) slope leading to apex, (neg) slope following apex, Tpeak-Send, rrms(Upeak-Send), width of Uwave , Upeak-Send]
depSlope = 2;
pinit(1) = depSlope;                                            % initial slope estimation depolarization
pinit(2) = 0.01;                                                % the initial value of the Tdom like function; = approx: y(0)
pinit(3) = -0.03;                                               % determines the (positive) slope leading up to the apex
pinit(4) = 0.06;                                                % determines the (negative) slope following the apex
pinit(5) = SPECS.time_apexT - SPECS.time_Jpoint;                % time between peak T-wave and end QRS
pinit(6) = SPECS.rrms(SPECS.time_apexU - SPECS.time_Jpoint);    % rrms between peak U-wave and end QRS
pinit(7) = 50;                                                  % width (sigma) of the U wave
pinit(8) = SPECS.time_apexU - SPECS.time_Jpoint;                % time between peak U-wave and end QRS
if funtype == 6, pinit = pinit(1:5); end                        % if funtype = 6, forget U-wave! 

tdom = lowpassma(SPECS.tdom(SPECS.time_Jpoint:SPECS.endtwave-SPECS.onsetqrs),20);   % lowpass filter rms for (Jpoint till end signal) minus timepoint onset QRS
tdom = baselinecor(tdom,find(tdom==min(tdom(1:200))),length(tdom))';                % baselinecorrect tdom from minimum in first 200 samples tdom till end
x    = (1:length(tdom))';

% compute parameters for TMP from fit to dominant T wave signal from which tdom is constructed:
% fit analytical function (product of logistics) to tdomS
paramsCumsum    = quamin_pvd(x,tdom,pinit,funtype);                             % input: [size rms, rms, initial estimates, funtype]: estimate measured tdom parameters
paramsCumsum(1) = max(depSlope,paramsCumsum(1));                                % if paramsCumsum(1) is lower than initial depSlope, use depSlope!
tdomFit         = rfunc(x, paramsCumsum, 0, funtype);                           % construct Tdomfit with estimated parameters
tdomFit         = max(SPECS.rrms(SPECS.time_Jpoint:end))/max(tdomFit)*tdomFit;  % scale-factor.

% construct figure for dominant t-wave [actual and fitted]
h302 = figure(302);
set(h302,'Name','Dominant T');
clf
plot(x,tdom,'b')
hold on
plot(x,tdomFit,'k')

% compute parameters for TMP from fit to derived TMP waveform
APtwave     = 1 - cumsum(tdom);                     % 1 minus cumulative [rms of signal BSM]
APtwave     = APtwave - APtwave(end);               % set last value to baseline [0]
APtwave     = APtwave / max(APtwave);               % normalize
APmean      = 1 - cumsum(tdomFit);                  % 1 minus cumulative [estimation signal BSM]
APmean      = APmean - APmean(end);                 % set last value to baseline [0]
APmean      = APmean / max(APmean);                 % normalize
pinit(1)    = max(depSlope,paramsCumsum(1));        % see line 311??
pinit(2)    = paramsCumsum(2);                      % change pinit(2) to final value
pinit(3)    = paramsCumsum(3);                      % change pinit(3) to final value
pinit(4)    = paramsCumsum(4);                      % change pinit(4) to final value

if ~SPECS.useCumsum, % if SPECS.useCumsum = 1, do not perform this part!
    % fit analytical function (logistics) to tdom
    parms       = quamin_pvd(x,APtwave,pinit,funtype);  % finetune estimated parameters with cumulative signal
    APfit       = rfunc(x,parms,0,funtype);
    APfit       = APfit - APfit(end);                   % set last value to baseline [0]
    APfit       = APfit / max(APfit);                   % normalize
    parms(5)    = parms(5) - (SPECS.time_apexT - SPECS.time_Jpoint);
end

paramsCumsum(5) = paramsCumsum(5) - (SPECS.time_apexT - SPECS.time_Jpoint);

if SPECS.useCumsum, p = paramsCumsum; else p = parms; end

% fill in final estimated values, based on Twave fit:
SPECS.depSlope      = p(1);
SPECS.initialSlope  = p(2);
SPECS.plateauslope  = p(3);
SPECS.repslope      = p(4);
SPECS.repCorrection = p(5);

% create examples of activation for visualization
dep                 = 10;
rep                 = SPECS.time_apexT + dep;
t                   = 1:length(tdom)+100;
S                   = getSmode(t, dep,    rep,       SPECS,4);
S1                  = getSmode(t, dep,    rep - 100, SPECS,4);
S2                  = getSmode(t, dep+100,rep - 100, SPECS,4);

% construct figure:
figure(303);
clf;

subplot(2,1,1)
plot(APtwave,'b','linewidth',1)
hold on;
plot(APmean,'r','linewidth',1);
axis tight
legend('1-int(Tdom)','1-int(fitted Tdom)')

subplot(2,1,2)
plot(t,S(1,:),'r','linewidth',1)
hold on
plot(t,S1(1,:),'r','linewidth',1)
plot(t,S2(1,:),'r','linewidth',1)
axis tight