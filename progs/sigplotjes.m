% sigplotjes.m
% dedicated version of sigplot for illustrations in papers
% 20050123
% function sigplotjes(PHI,nrows,ncols,LABS,funscal,col,tbeg,tend,linw)
function sigplotjes(PHI,nrows,ncols,LABS,funscal,col,tbeg,tend,linw)
if nargin<0, linw=1; end
[nlds nt]=size(PHI);
LAY=[];
nsig=1;
LAY(1,1:3)=[nrows ncols 0];
for i=2:nrows+1,
    for j=1:ncols,
       nsig=nsig+1;
       LAY(nsig,1:3)=[i-1 j nsig-1];
   end
end


nplts=nlds;
if ~exist('label'), label=' '; end
if ~exist('funscal'), funscal=1; end
if ~exist('col'), col='blue'; end
if ~exist('leadnos'), leadnos=1; end
if ~exist('markt'), markt=0; end

grid=LAY(1,:);
axis([0 grid(2) 0 2*grid(1)] )
axis manual
axis('off')
hold on
t=1:nt;
t=t*.8/nt;

for i=2:nplts+1,
   j=LAY(i,3);
   yshift(j)=2*grid(1)+1-2*LAY(i,1);
   xshift(j)=LAY(i,2)-1;
   plot(t(tbeg:tend)+xshift(j),funscal*PHI(j,tbeg:tend)+yshift(j),col,'linewidth',linw)
   plot([t(1)+xshift(j) t(nt)+xshift(j)],[yshift(j) yshift(j)],':k')
   text(xshift(j)+.7,yshift(j)+.5,LABS(j,:),'fontsize',12);
end

extr=extremes(PHI(:,tbeg:tend));
sprintf('Extremes: min %5.3f  at: %5d %5d; max %5.3f  at: %5d  %5d', extr)

next=0;
if next==1,
    plot(t(extr(3))+xshift(extr(2)),funscal*extr(1)+yshift(extr(2)),'*b')
    plot(t(extr(6))+xshift(extr(5)),funscal*extr(4)+yshift(extr(5)),'*r')
end

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
zm=0;
if zm==0,
    funlabel=uicontrol('style','text','units','norm','position',[.6 .95 .4 .05],...
   'string',label,'fontsize',8);
else,
    funlabel=uicontrol('style','text','units','norm','position',[.6 .95 .4 .05],...
   'string',['zm;' label],'fontsize',8);
end

xl=[0   0  0.025  0];
yl=[1   2  1.8   1.8];
x2=[0.1 0.1     0    0.8 0.7 0.7];
y2=[0   0.025   0  0 0.025 0];
set(line(xl,yl),'color','k')

set(line(x2,y2),'color','k')

%text(0.05, 1.9,sprintf('%0.3f mV',1/funscal))
text(0.05, 1.9,'100')

text(0.2, 0.15,sprintf('%d %s',nt,'ms'))

