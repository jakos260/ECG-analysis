function v12holter(varargin) %(t,PHI)

% plot standard leads

VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols=['b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];

% cols=['k  ';'k: ';'k--';'k-.';];

% sampT  is expected in hours
warning off
paperspeed=2.5;
do9=0;
dowct =0;
tmax=0;
marks=0;
info=[];
maxphi=0;
minphi=0;
sampT=1/60;
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
tmax = ceil(tmax*sampT);

if length(maxphi)>1
    minphi=maxphi(1);
    maxphi=maxphi(2);
elseif maxphi~=0
    minphi=-maxphi;
end

%% wct reference and determine maxphi
wct = 1:3;
vs =  [wct 4:9];

for i=1:nmap
    A=PHI{i};
    if dowct
        A = A - ones(size(A,1),1)*mean(A(wct,:));
    end
    PHI{i} = A;
    if maxphi==0
        maxphi=max(maxphi,max(PHI{i}(:)));
        minphi=min(minphi,min(PHI{i}(:)));
    end
end

maxphi = max(abs([maxphi minphi]));
maxphi = ceil(maxphi * 2)/2;
minphi = -maxphi;

% if paperspeed<=50
%     tmax=ceil(tmax/paperspeed)*paperspeed;
% end



% clf
linew=2;
fs=9;

fweight='demi';
clf

if do9
    ny=3;nx=3;
    x=0;
else
    x=0;
    ny=3;nx=4;
    textI=['I  ';'II ';'III'];
    for k=1:3
        y=3-icyc(k,3);
        axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
        for i=nmap:-1:1     
            hold on
            if size(PHI{i},1) == 12
                lead1 = PHI{i}(1,:);
                lead2 = PHI{i}(2,:);
                lead3 = PHI{i}(3,:);
            else
                phi1=PHI{i}(wct(1),:);
                phi2=PHI{i}(wct(2),:);
                phi3=PHI{i}(wct(3),:);
                lead1 = phi2-phi1;
                lead2 = phi3-phi1;
                lead3 = phi3-phi2;
            end
            color = cols(i,:);
            if ~isempty(mycolor)
                color = mycolor;
            end
            if ~bw
                if k==1
                    plot(sampT*(0:length(lead1)-1), lead1/1.5,color,'Linewidth',linew);
                elseif k==2
                    plot(sampT*(0:length(lead2)-1), lead2/1.5,color,'Linewidth',linew);
                else
                    plot(sampT*(0:length(lead3)-1), lead3/1.5,color,'Linewidth',linew);
                end
            else
                if k==1
                    plot(sampT*(0:length(lead1)-1), lead1,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                elseif k==2
                    plot(sampT*(0:length(lead2)-1), lead2,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                else
                    plot(sampT*(0:length(lead3)-1), lead3,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                end
            end
        end
        
        axis([0,tmax,minphi,maxphi]); axis off
        text(tmax*0.65,0.9*maxphi,textI(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
        ecgraster('Amplification',Amplification,'SampleRate',60*sampT,'Paperspeed',paperspeed)
    end
    x=1;
end

for k=1:9
    y=3-icyc(k,3);
    axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
    if mod(k,ny)==0, x=x+1; end
    for i=nmap:-1:1
        doCorrect =1;
        if size(PHI{i},1) == 12
            doCorrect  = 0;
            PHI9 = PHI{i}(4:end,:);
        else
            PHI9 = PHI{i};
        end
        V=vs(k);
        phi=PHI9(V,:);
        if k >= 1 && k <=3 && doCorrect
            phi= phi* 1.5;
        end
        color = cols(i,:);
        if ~isempty(mycolor)
            color = mycolor;
        end
        if ~bw
            plot(sampT*(0:length(phi)-1), phi,color,'Linewidth',linew);
        else
            plot(sampT*(0:length(phi)-1), phi,styles(icyc(i,length(styles)),:),'Linewidth',linew);
        end
        hold on;
    end
    
    axis([0,tmax,minphi,maxphi]); axis off
    text(tmax*0.5,0.9*maxphi,VLab(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
    ecgraster('Amplification',Amplification,'SampleRate',60*sampT,'Paperspeed',paperspeed)
    
end
annotation('textbox',[0.65 0.01 0.6 0.025],'string',[num2str(paperspeed) ' mm/h  ' num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');
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
