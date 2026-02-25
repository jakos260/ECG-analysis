function leadv12(varargin) %(t,PHI)

% plot standard leads
KI_9=[ -0.5267    0.1634   -0.2867   -0.1300    0.0500   -0.0100    0.1400    0.0600    0.5400;   % Vx
       -0.8633   -0.0733    0.9266    0.0600   -0.0200   -0.0500    0.0600   -0.1700    0.1300;   % Vy
        0.3300    0.3200   -0.0200   -0.4300   -0.0600   -0.1400   -0.2000   -0.1100    0.3100;];  % Vz
    
    
VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols=['b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];
warning off
paperspeed=50;
do9=0;
dowct =0;
tmax=0;
marks=0;
info=[];
maxphi=0;
minphi=0;
sampT=1/1000;
mycolor=[];
bw=0;
nmap=0;
Amplification=1;
PHI=[];
if length(varargin) < 1
    error('This routine needs at least two parameters');
else
    pp=1;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case '9leads'
                    do9=varargin{pp+1};pp=pp+2;
                case 'marks'
                    marks=varargin{pp+1};pp=pp+2;
                case 'do9'
                    do9=varargin{pp+1};pp=pp+2;
                case 'dowct'
                    dowct=varargin{pp+1};pp=pp+2;
                case 'info'
                    info=varargin{pp+1};pp=pp+2;
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
                case 'color'
                    mycolor = varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        else
            if iscell(varargin{pp})
                
                for i=1:length(varargin{pp})
                    nmap = nmap+1;                 
                    PHI{nmap} = varargin{pp}{i};
                    tmax = max(tmax,size(PHI{pp},2)-1);                    
                end
            else
                nmap = nmap+1;
                
                PHI{nmap} = varargin{pp};
                tmax = max(tmax,size(PHI{pp},2)-1);
                
            end
            
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
wct = 1:3;

for i=1:nmap
    if size(PHI{i},1) == 9
        VCG{i} = KI_9 * PHI{i};
    elseif size(PHI{i},1) == 12
        VCG{i} = KI_9 * PHI{i}(4:end,:);
    else
        error('not 12 leads');
    end
    maxphi=max(maxphi,max(VCG{i}(:)));
    minphi=min(minphi,min(VCG{i}(:)));
end

maxphi = max(abs([maxphi minphi]));
maxphi = ceil(maxphi * 2)/2;
minphi = -maxphi;




% clf
linew=1;
fs=9;

fweight='demi';
clf


ny=1;nx=3;y=0;
textI=['X';'Y';'Z'];
for k=1:3  
    x=k-1;
    axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
    for i=nmap:-1:1
        hold on
        lead1 = VCG{i}(1,:);
        lead2 = VCG{i}(2,:);
        lead3 = VCG{i}(3,:);
        color = cols(i,:);
        if ~isempty(mycolor)
            color = mycolor;
        end
        if k==1
            plot(lead2, lead1,color,'Linewidth',linew);
        elseif k==2
            plot(lead2, lead3,color,'Linewidth',linew);
        else
            plot(lead1, lead3,color,'Linewidth',linew);
        end
    end

    axis([minphi,maxphi,minphi,maxphi]); axis off
    text(tmax*0.65,0.9*maxphi,textI(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
    ecgraster('SampleRate',1000/sampT,'Paperspeed',paperspeed,'Amplification',Amplification)
end



annotation('textbox',[0.65 0.01 0.6 0.025],'string',[ num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');
%%
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
    xw=diff(xlim)/Amplification*Scale;
    yh=diff(ylim)/Amplification*Scale;
    pos=get(gca,'position');
    set(gca,'xlim',xl,'ylim',yl,'Position',[pos(1),pos(2),xw,yh]);
    set(gca,'Units',oldunits);
end
daspect([0.5*Amplification,0.5*Amplification,1]);

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
if diff(xl)/(.1*Amplification)<1000
    for k=xl(1):.1*Amplification:xl(2);
        plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',.5)
    end
    %    for k=xl(1):5/Paperspeed:xl(2);
    %       plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',1)
    %    end
    for k=xl(1):.5*Amplification:xl(2);
        plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',1.5)
    end
    
else
    disp(['Too many lines in ',mfilename,', possibly wrong sample freq']);
end
set(UseAxis,'children',[c0;setdiff(get(UseAxis,'children'),c0)]);
