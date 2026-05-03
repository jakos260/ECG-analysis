% sigplot voor matlab
% function sigplot(PHI,label,LAY,funscal,col,leadnos,markt,listextr,linew)
% 20060407

% function extr=sigplot(PHI,label,LAY,funscal,col,leadnos,markt,listextr,linew)

function extr=sigplot(varargin)
	
global sinter;



if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	PHI=varargin{1};
	label=' '; 
	funscal=1;
	col='blue';
	leadnos=0; 
	markt=0; 
	listextr=0;
	sampT=1;
	linew=1;
	
	pp=2;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'label'
					label=varargin{pp+1};pp=pp+2;
				case 'funscal'
					funscal=varargin{pp+1};pp=pp+2;
				case 'lay'
					LAY=varargin{pp+1};pp=pp+2;
				case 'sampt'
					sampT==varargin{pp+1};pp=pp+2;
				otherwise
					switch pp
						case 2 
							label=varargin{pp};pp=pp+1;
						case 5
							col = varargin{pp};pp=pp+1;
						otherwise
							error('unknown parameter');
					end
			end
		else
			switch pp
				case 3
					LAY= varargin{pp};pp=pp+1;
				case 4
					funscal= varargin{pp};pp=pp+1;
				case 5
					col= varargin{pp};pp=pp+1;
				case 6 
					leadnos = varargin{pp};pp=pp+1;
				case 7 
					markt= varargin{pp};pp=pp+1;
				case 8
					listextr= varargin{pp};pp=pp+1;
				case 9
					linew= varargin{pp};pp=pp+1;
				otherwise
					error('unknown parameter');
			end
		end
	end
end




[nlds nt]=size(PHI);
[nplts jdum]=size(LAY); % note: nplts may be smaller than nlds
nplts=nplts-1;

if nplts==12,
   % test for augmented leads
   augm=0;
   if norm(PHI(5,:)-PHI(4,:))/norm(PHI(1,:))>1.1,
      augm=1;
      LAB12=['  I'; ' II'; 'III'; 'AVR'; 'AVL'; 'AVF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
   else,
      LAB12=['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
   end
end
if size(LAY,1)==9,
    LAB12=['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; 'V2s'; 'V1s'; 'Vrs'; 'V2R';];
end
 LAB9=[' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];

% test for zero mean
zm=0;
allmean=mean(mean(PHI));
if abs(allmean) < 0.002,
   zm=1;
end

grid=LAY(1,:);
axis([0 grid(1) 0 2*grid(2)] )
axis manual
axis('off')
hold on
t=1:sampT:nt*sampT;
ngridcols=LAY(1,1);
t=t*.8/nt*ngridcols/(ngridcols-0.2);

for i=2:nplts+1
   j=LAY(i,3);
   yshift(j)=2*grid(2)+1-2*LAY(i,2);
   xshift(j)=LAY(i,1)-1;
   plot(t+xshift(j),funscal*PHI(j,:)+yshift(j),col,'linewidth',linew)
   plot([t(1)+xshift(j) t(nt)+xshift(j)],[yshift(j) yshift(j)],':k')
   if leadnos==1,
       if nlds==12,
           text(xshift(j)+.4,yshift(j)+.25,LAB12(j,:));
       end
       if nlds==9,
           text(xshift(j)+.4,yshift(j)+.25,LAB9(j,:));
       else,
           if nlds~=12,
          text(xshift(j)+.4,yshift(j)+.25,num2str(j)); 
      end
       end
   end
end

extr=extremes(PHI);
% disp(sprintf('Extremes: min %5.3f  at: %5d %5d; max %5.3f  at: %5d  %5d', extr));
plot(t(extr(3))+xshift(extr(2)),funscal*extr(1)+yshift(extr(2)),'*b')
plot(t(extr(6))+xshift(extr(5)),funscal*extr(4)+yshift(extr(5)),'*r')

if markt>0 & markt<=nt,
   tmark=t(markt);
   MARKS=[];
   for i=2:nplts+1,
      j=LAY(i,3);
      ymark=funscal*PHI(j,markt)+yshift(j);
      YMARK=[ymark-.1; ymark+.1];
      TMARK=tmark;
      TMARK=ones(2,1)*TMARK;
      MARKS=[MARKS plot(TMARK+xshift(j),YMARK,'r')];
    end
end

if zm==0,
    funlabel=uicontrol('style','text','units','norm','position',[.6 .95 .4 .05],...
   'string',label,'fontsize',8);
else,
    funlabel=uicontrol('style','text','units','norm','position',[.6 .95 .4 .05],...
   'string',['zm;' label],'fontsize',8);
end

xl=[0 0  1];
yl=[1 0 0 ];
set(line(xl,yl),'color','k')
text(0.05, 0.6,sprintf('%0.3f mV',1/funscal))

if sinter==[], sinter=1;  end
text(-LAY(1,1)/10, -.25,sprintf('%d %s %0.3f %s',nt,'samples; ',(nt-1)*sinter/1000,'s'))

if listextr~=0,
   laag=min(min(PHI));
   hoog=max(max(PHI));
   text(grid(1)/2, -.7,sprintf('%s %0.2f  to  %0.2f mV',' range ',laag,hoog))
end
