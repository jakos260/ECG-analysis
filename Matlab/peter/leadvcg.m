function leadvcg(varargin) %(t,PHI)

% plot standard leads

VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols=['b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];

% cols=['k  ';'k: ';'k--';'k-.';];


warning off
paperspeed=50;
tmax=0;
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
                case 'max'
                    maxphi=varargin{pp+1};pp=pp+2;
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

elec12ECG=[11.4824 -200.0000  174.4150;...
          -4.9266  199.2190  182.4210;...
          40.9138   47.1829 -285.5690;...
          96.6231  -24.9980   38.4330;...
          91.1554   30.4652   38.4342;...
          92.7157   83.5986   -1.5735;...
          73.9623  123.4420    6.4327;...
          28.6549  146.8750   -1.5688;...
         -33.0585  176.5600   -1.5703];
pos = [39.9854   21.9069    5.5178];

%% wct reference and determine maxphi
wct = 1:3;
vs =  [wct 4:9];

for i=1:nmap
    A=PHI{i};
    if size(A,1) == 9
        A = [A(1:3,:)/1.5 ; A(4:end,:)];
    elseif size(A,1) == 12
        A = [A(4:6,:)/1.5 ; A(7:end,:)];
    else
        A = [A(4:6,:)/1.5 ; A(7:end,:)];
    end
        
    [vcg,korsvcg,dowervcg]=ecgvector(A,elec12ECG,pos);
    disp(['kors area ' num2str(norm3d(sum(korsvcg))) '   dower area ' num2str(norm3d(sum(dowervcg)))])
    PHI{i} = korsvcg;
    if maxphi==0
        maxphi=max(maxphi,max(PHI{i}(:)));
        minphi=min(minphi,min(PHI{i}(:)));
    end
end

maxphi = max(abs([maxphi minphi]));
maxphi = ceil(maxphi * 2)/2;
minphi = -maxphi;

if paperspeed<=50
    tmax=ceil(tmax/paperspeed)*paperspeed;
end



% clf
linew=1;
fs=9;

fweight='demi';
clf

x=0;
ny=2;nx=3;
textI=['X';'Y';'Z'];
for x=0:2
    y=1;
    axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
    for i=nmap:-1:1     
        hold on
        vcg = PHI{i};
        color = cols(i,:);
        if ~isempty(mycolor)
            color = mycolor(i,:);
        end
        plot(sampT*(0:size(vcg,1)-1), vcg(:,x+1),color,'Linewidth',linew);
    end
    axis([0,tmax,minphi,maxphi]); axis off
    title(textI(x+1,:))
%     text(tmax*0.65,0.9*maxphi,textI(x+1,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
    ecgraster('SampleRate',1000/sampT,'Paperspeed',paperspeed,'Amplification',Amplification)
end

y=0;
for x=0:2
    axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
    vcgraster('Scale',maxphi)
    for i=nmap:-1:1     
        hold on
        vcg = PHI{i};
        color = cols(i,:);
        if ~isempty(mycolor)
            color = mycolor(i,:);
        end
         if x==0
            plot(vcg(:,2), -vcg(:,1),color,'Linewidth',linew);
            title('horizontal');
            axis off equal
        elseif x==1
            plot(vcg(:,1), vcg(:,3),color,'Linewidth',linew);
            title('right sagital');
            axis off equal
        elseif x==2
            plot(vcg(:,2), vcg(:,3),color,'Linewidth',linew);
            title('frontal');
            axis off equal
        end
    end


end

annotation('textbox',[0.65 0.01 0.6 0.025],'string',[num2str(paperspeed) ' mm/s  ' num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');
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

%%
function vcgraster(varargin)
Amplification=1; % 1mV/mm
Scale=1;
UseAxis=gca;
ProcessVarargin(varargin);

axes(UseAxis);
c0=get(UseAxis,'children');
hold on
% daspect([0.5*Amplification,0.5*Amplification,1]);
axis([-Scale,Scale,-Scale,Scale])
yl=ylim;
xl=xlim;
if diff(yl)/(.1*Amplification)<1000
    for k=yl(1):.1*Amplification:yl(2)
        plot(xl,k*[1,1],'color',[1,.5,.5],'linewidth',.5)
        plot(k*[1,1],xl,'color',[1,.5,.5],'linewidth',.5)
    end
    for k=yl(1):.5*Amplification:yl(2)
        plot(xl,k*[1,1],'color',[1,.5,.5],'linewidth',1.5)
        plot(k*[1,1],xl,'color',[1,.5,.5],'linewidth',1.5)
    end
else
    disp(['Too many lines in ',mfilename,', possibly wrong amplitude']);
end

set(UseAxis,'children',[c0;setdiff(get(UseAxis,'children'),c0)]);

