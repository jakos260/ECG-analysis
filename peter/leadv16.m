function leadv16(varargin) %(t,PHI)

% plot standard leads from the nijmegen leads system

VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols=['b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];

paperspeed=50;
do9=0;
dowct = 0;
tmax=0;
marks=0;
info=[];
maxphi=0;
minphi=0;
sampT=1/1000;
bw=0;
fixedcolor =[];
Amplification=1;
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
				case 'leadsys'
					leadsys=varargin{pp+1};pp=pp+2;					
				case 'sampt'
					sampT=varargin{pp+1};pp=pp+2;					
				case 'paperspeed'
					paperspeed=varargin{pp+1};pp=pp+2;					
				case 'amplification'
					Amplification=varargin{pp+1};pp=pp+2;							
                case 'color'
                    fixedcolor = varargin{pp+1};pp=pp+2;							
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
tmax=tmax*sampT+10*sampT;

if length(maxphi)>1
	minphi=maxphi(1);
	maxphi=maxphi(2);
elseif maxphi~=0
	minphi=-maxphi;
end

if ~isempty(strfind(leadsys,'nim')) || ~isempty(strfind(leadsys,'nijmegen'))
	wct=[1 2 65];
	vs=[wct 19 26 34 41 48 54];	
	if size(PHI_1,1)==64
		for i=1:nmap
			eval(['A=PHI_' num2str(i) ';']); 
			A=zeromean([A;-(A(1,:)+A(2,:))]);
			eval(['PHI_' num2str(i) '=A;']);
		end
    end
    dowct = 1;
elseif  strfind(leadsys,'ams')
	wct=[64 63 65];
	vs=[wct 12 18 25 31 40 45];
    dowct = 1;
elseif  strfind(leadsys,'arnhem')
	wct=[97 98 99];
	vs=[wct 20 28 37 45 53 61];
    dowct = 1;
elseif ~isempty(strfind(leadsys,'12lead')) || ~isempty(strfind(leadsys,'9lds'))
	wct=[1 2 3];
	vs=[wct 4 5 6 7 8 9];
elseif strfind(leadsys,'pig07')
	wct=[31 63 64];
	vs=[wct 21 27 36 42 46 49];   
    dowct = 1;
end
%% wct reference and determine maxphi

for i=1:nmap
	eval(['A=PHI_' num2str(i) ';']); 
    if dowct
        A=A-ones(size(A,1),1)*mean(A(wct,:));
    end
	eval(['PHI_' num2str(i) '=A;']);
	if maxphi==0
		eval(['maxphi=max(max(max(maxphi,((PHI_' num2str(i) '(vs,:))))));']);
		eval(['minphi=min(min(min(minphi,((PHI_' num2str(i) '(vs,:))))));']);
	end
end

% maxphi=round(maxphi);
% minphi=round(minphi);
if paperspeed<=50
	tmax=ceil(tmax/paperspeed)*paperspeed;
end

if ~isempty(info)
	if marks>0
		eval(['PHI_m=PHI_' num2str(marks) ';']);	
	else
		marks=1;			
	end
elseif marks>0
	eval(['PHI_m=PHI_' num2str(marks) ';']);
	info=getAmplInfos(PHI_m);
end

% clf
pos=get(gcf,'Position');
maxsize=get(0,'ScreenSize');
% if do9
% 	set(gcf,'Position',[pos(1)  min(pos(2),520) 520   520]);
% else
% % 	set(gcf,'Position',[pos(1)  min(pos(2),maxsize(4)-1.33333*520-150) 1.33333*520   520]);
% % 	set(gcf,'Position',[pos(1)  pos(2) 770 940]);	
% 	set(gcf,'Position',[pos(1)  pos(2) 485 600]);		
% end
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
			eval(['PHI=PHI_' num2str(i) ';']); hold on
            if isempty( fixedcolor )
                myColor = cols(i);
            else
                myColor = fixedcolor;
            end
			phi1=PHI(wct(1),:);
			phi2=PHI(wct(2),:);
			phi3=PHI(wct(3),:);				
			if ~bw
				if k==1
					plot(sampT*(0:length(phi2)-1), phi2-phi1,myColor,'Linewidth',linew);
				elseif k==2
					plot(sampT*(0:length(phi3)-1), phi3-phi1,myColor,'Linewidth',linew);
				else
					plot(sampT*(0:length(phi2)-1), phi3-phi2,myColor,'Linewidth',linew);				
				end
			else
				if k==1
					plot(sampT*(0:length(phi2)-1), (phi2-phi1)/1.5,styles(icyc(i,length(styles)),:),'Linewidth',linew);
				elseif k==2
					plot(sampT*(0:length(phi3)-1), (phi3-phi1)/1.5,styles(icyc(i,length(styles)),:),'Linewidth',linew);
				else
					plot(sampT*(0:length(phi2)-1), (phi3-phi2)/1.5,styles(icyc(i,length(styles)),:),'Linewidth',linew);				
				end
			end
		end
		
		
		% make sure the first signal is always seen
		eval(['PHI=PHI_' num2str(1) ';']);	hold on; 
		phi1=PHI(wct(1),:);
		phi2=PHI(wct(2),:);
		phi3=PHI(wct(3),:);				
		axis([0,tmax,minphi,maxphi]); axis off
		text(tmax*0.65,0.9*maxphi,textI(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);	
		ecgraster('SampleRate',1000/sampT,'Paperspeed',paperspeed,'Amplification',Amplification)			
	end
	x=1;
end

for k=1:9
	y=3-icyc(k,3);
	axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
	if mod(k,ny)==0, x=x+1; end
	for i=nmap:-1:1
		eval(['PHI=PHI_' num2str(i) ';']); hold on
        if isempty( fixedcolor )
            myColor = cols(i);
        else
            myColor = fixedcolor;
        end

		V=vs(k);
		phi=PHI(V,:);
        if k >= 1 && k <=3
            phi= phi* 1.5;
        end
		if k==6 %V3
			if ~bw
				plot(sampT*(0:length(PHI(V,:))-1), mean(PHI(V-1:V,:)),myColor,'Linewidth',linew);
			else
				plot(sampT*(0:length(PHI(V,:))-1), mean(PHI(V-1:V,:)),styles(icyc(i,length(styles)),:),'Linewidth',linew);				
			end
			if marks~=0
				plot(sampT*mean(info(V-1:V,5)),mean(mean(PHI_m(V-1:V,info(V-1:V,5)))),'yo','MarkerFacecolor','y');
				plot(sampT*mean(info(V-1:V,6)),mean(mean(PHI_m(V-1:V,info(V-1:V,6)))),'ro','MarkerFacecolor','r');			
				plot(sampT*mean(info(V-1:V,7)),mean(mean(PHI_m(V-1:V,info(V-1:V,7)))),'go','MarkerFacecolor','g');			
				plot(sampT*mean(info(V-1:V,8)),mean(mean(PHI_m(V-1:V,info(V-1:V,8)))),'ko','MarkerFacecolor','k');			
			end
		elseif k<=3
			if ~bw
				plot(sampT*(0:length(phi)-1), phi,myColor,'Linewidth',linew);
			else
				plot(sampT*(0:length(phi)-1), phi,styles(icyc(i,length(styles)),:),'Linewidth',linew);
			end
			if marks ~=0
				plot(sampT*info(V,5),PHI_m(V,info(V,5)),'yo','MarkerFacecolor','y');
				plot(sampT*info(V,6),PHI_m(V,info(V,6)),'ro','MarkerFacecolor','r');			
				plot(sampT*info(V,7),PHI_m(V,info(V,7)),'go','MarkerFacecolor','g');			
				plot(sampT*info(V,8),PHI_m(V,info(V,8)),'ko','MarkerFacecolor','k');			
			end
		else
			if ~bw
				plot(sampT*(0:length(PHI(V,:))-1), PHI(V,:),myColor,'Linewidth',linew);
			else
				plot(sampT*(0:length(PHI(V,:))-1), PHI(V,:),styles(icyc(i,length(styles)),:),'Linewidth',linew);
			end
			if marks ~=0
				plot(sampT*info(V,5),PHI_m(V,info(V,5)),'yo','MarkerFacecolor','y');
				plot(sampT*info(V,6),PHI_m(V,info(V,6)),'ro','MarkerFacecolor','r');			
				plot(sampT*info(V,7),PHI_m(V,info(V,7)),'go','MarkerFacecolor','g');			
				plot(sampT*info(V,8),PHI_m(V,info(V,8)),'ko','MarkerFacecolor','k');			
			end
		end
	end
	% make sure the first signal is always seen
	eval(['PHI=PHI_' num2str(1) ';']);	hold on; 
	axis([0,tmax,minphi,maxphi]); axis off
	text(tmax*0.5,0.9*maxphi,VLab(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
	ecgraster('Amplification',Amplification,'SampleRate',1000/sampT,'Paperspeed',paperspeed)	
end
% pos=get(gcf,'Position');
% set(gcf,'Position',[pos(1),pos(2) 400 840])
annotation('textbox',[0.65 0.01 0.6 0.025],'string',[num2str(paperspeed) ' mm/s  ' num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');
% if tmax <500
% 	AxisLegend(t,[-maxphi maxphi],100,1);
% else
% 	AxisLegend(t,[-maxphi maxphi],200,1);
% end

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
