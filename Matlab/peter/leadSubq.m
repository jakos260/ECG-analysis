function leadSubq(varargin) %(t,PHI)

VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols=['b  ';'k  ';'r  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];

paperspeed=50;
do9=0;
dowct =0;
tmax=0;
marks=0;
info=[];
maxphi=0;
minphi=0;
sampT=1/1000;
bw=0;
Amplification=1;
if length(varargin) < 1
    error('This routine needs at least two parameters');
else
    pp=1;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'max'
                    maxphi=varargin{pp+1};pp=pp+2;
                case 'bw'
                    bw=varargin{pp+1};pp=pp+2;
                case 'sampt'
                    sampT=varargin{pp+1};pp=pp+2;
                case 'paperspeed'
                    paperspeed=varargin{pp+1};pp=pp+2;
                case 'amplification'
                    Amplification=varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        else
            eval(['PHI_' num2str(pp) '=varargin{' num2str(pp) '};']);
            eval(['tmax=max(size(PHI_' num2str(pp) ',2)-1,tmax);']);
            nmap=pp;
            pp=pp+1;
        end
    end
end
sampT=1000*sampT;
tmax=tmax*sampT+1*sampT;

if length(maxphi)>1
    minphi=maxphi(1);
    maxphi=maxphi(2);
elseif maxphi~=0
    minphi=-maxphi;
end

%% wct reference and determine maxphi


for i=1:nmap
    eval(['A=PHI_' num2str(i) ';']);
    if dowct
        A = A - ones(size(A,1),1)*mean(A(wct,:));
    end
    eval(['PHI_' num2str(i) '=A;']);
    if maxphi==0
        eval(['maxphi=max(max(max(maxphi,((PHI_' num2str(i) '(:,:))))));']);
        eval(['minphi=min(min(min(minphi,((PHI_' num2str(i) '(:,:))))));']);
    end
end

% maxphi=round(maxphi);
% minphi=round(minphi);
if paperspeed<=50
    tmax=ceil(tmax/paperspeed)*paperspeed;
end


linew=1.5;
fs=9;

fweight='demi';
clf
ny = 1;
nx = size(A,1);
for k=1:nx
    x = k-1;
    axes('Position',[x/nx+0.005,0,1/nx-0.015,0.99*ny]);
    for i=nmap:-1:1
        eval(['PHI=PHI_' num2str(i) ';']); hold on
        if ~bw
            plot(sampT*(0:length(PHI(k,:))-1), PHI(k,:),cols(i,:),'Linewidth',linew);
        else
            plot(sampT*(0:length(PHI(k,:))-1), PHI(k,:),styles(icyc(i,length(styles)),:),'Linewidth',linew);
        end
    end
    % make sure the first signal is always seen
    eval(['PHI=PHI_' num2str(1) ';']);	hold on;
    axis([0,tmax,minphi,maxphi]); axis off
    text(tmax*0.5,0.9*maxphi,['lead ' num2str(k)],'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
    hold on;
    ecgraster('Amplification',Amplification,'SampleRate',1000/sampT,'Paperspeed',paperspeed)
end
annotation('textbox',[0.65 0.03 0.6 0.025],'string',[num2str(paperspeed) ' mm/s  ' num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');


function ecgraster(varargin)
Amplification=1; % 1mV/mm
Paperspeed=25; % 25 mm/s
Calibrate='off';
Scale=1;
UseAxis=gca;
SampleRate=1000; %
ProcessVarargin(varargin);
Paperspeed=Paperspeed/SampleRate;

axes(UseAxis);
c0=get(UseAxis,'children');
hold on
if strcmpi(Calibrate,'on')
    oldunits=get(gca,'Units');
    set(gca,'Units','centimeters');
    yl=ylim;xl=xlim;
    xw=diff(xlim)*Paperspeed/10*Scale;
    yh=diff(ylim)/Amplification*Scale;
    pos=get(gca,'position');
    set(gca,'xlim',xl,'ylim',yl,'Position',[pos(1),pos(2),xw,yh]);
    set(gca,'Units',oldunits);
end
daspect([5/Paperspeed,0.5*Amplification,1]);

yl=ylim;
xl=xlim;
if diff(yl)/(.1*Amplification)<1000
    for k=yl(1):.1*Amplification:yl(2);
        plot(xl,k*[1,1],'color',[1,.5,.5],'linewidth',.5)
    end
    for k=yl(1):.5*Amplification:yl(2);
        plot(xl,k*[1,1],'color',[1,.5,.5],'linewidth',1.5)
    end
else
    disp(['Too many lines in ',mfilename,', possibly wrong amplitude']);
end
if diff(xl)*Paperspeed<1000
    for k=xl(1):1/Paperspeed:xl(2);
        plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',.5)
    end
    %    for k=xl(1):5/Paperspeed:xl(2);
    %       plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',1)
    %    end
    for k=xl(1):5/Paperspeed:xl(2);
        plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',1.5)
    end
    
else
    disp(['Too many lines in ',mfilename,', possibly wrong sample freq']);
end
set(UseAxis,'children',[c0;setdiff(get(UseAxis,'children'),c0)]);
