% sigplot voor matlab
% function sigplot(PHI,label,LAY,funscal,col,signodes,markt,listextr,linew)
% 20170111
% if signodes==0         : no signal labels plotted
% if signodes==1         : signal labels plotted
% if signodes is a vector: its nsig  lements are used as labels

%

%function extr=sigplot(PHI,label,LAY,funscal,col,signodes,markt,listextr,linew)
function sigplot(PHI,label,LAY,funscal,col,signodes,markt,listextr,linew)
global sint;

grid=LAY(1,:);
LAY(1,:)=[];
nplts=size(LAY,1);

[nsigs nt]=size(PHI);
nplts=min(nplts,nsigs); % note: nplts may be smaller than nlds

if nargin<9, linew=1; end
if ~exist('label'), label=' '; end
if ~exist('funscal'), funscal=1; end
if ~exist('col'), col='blue'; end
if ~exist('signodes'), signodes=1; end
if ~exist('markt'), markt=0; end
if ~exist('listextr'), listextr=0; end

if max(size(signodes))==1,
    if signodes~=0,
        if nplts==12,
            % standard leads
            % test for augmented extremities
            augm=0;
            norm((PHI(5,:)-PHI(4,:)-PHI(1,:))/nt);
            
           
            
            
            if norm((PHI(5,:)-PHI(4,:)-PHI(1,:))/nt)>1.e-4,
                augm=1;
                LAB=['  I'; ' II'; 'III'; 'AVR'; 'AVL'; 'AVF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
            else,
                LAB=['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
            end
        end
        
        if size(LAY,1)==9,
            LAB=[' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
        end
        
        if nplts==13,
            % common reference leads; OACG
            LAB=[' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';' A1';' A2';' A3';' A4'; ];
        end
        %
        if nplts==16,
            % standard 12 signals plus four OACG leads
            LAB=['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';' A1';' A2';' A3';' A4'; ];
        end
        
        %LAB=[' 1';' 2'; ' 3'; ' 4'; ' 5'; ' 6'; ' 7'; ' 8'; ' 9' '10'; '11'; '12'; '13'; '14'; '15'];
    end
end

% zm=0
% allmean=mean(mean(PHI));
% if abs(allmean) < 0.002,
%     zm=1;
% end

axis([0 grid(1) 0 2*(grid(2)-1)] );
axis manual
axis('off')
hold on
t=1:nt;
ngridcols=grid(1);

t=t*.8/nt*ngridcols/(ngridcols-0.2);

for i=1:nplts,
    yshift(i)=2*grid(2)+1-2*LAY(i,2);
    xshift(i)=LAY(i,1)-1;
end

coll='b';

for i=1:nplts,
    isig=LAY(i,3);
    j=i;
    coll=col;  

    if ismember(isig,[1 2 19 26 33 41 48 54])   && nplts==64, % NIMDATA
        
    %if ismember(isig,[12 18 25 31 40 45 63 64 65]) && nplts>16, % nplts==65,   % ALBUM
        LAY(i,3);
        coll='r';
    end

    plot(t+xshift(j),funscal*PHI(isig,:)+yshift(j),coll,'linewidth',linew)
    
    if markt~=0,
        plot(t(markt)+xshift(j),funscal*PHI(isig,:)+yshift(j),'r','linewidth',linew)
    end
        
    plot([t(1)+xshift(j) t(nt)+xshift(j)],[yshift(j) yshift(j)],':k')
    
    if signodes==1
        if nplts<20 && nplts>1
            if ismember(size(LAB,1),[9 12 13 16])
                text(xshift(j)+.7,yshift(j)+.25,LAB(i,:));
            end
        else
            text(xshift(j)+.7,yshift(j)+.25,num2str(isig)); % for nim65
            %text(xshift(j)+.7,yshift(j)+.25,num2str(j));   % else
        end
    end
    
    if length(signodes)>1
        text(xshift(j)+.7,yshift(j)+.25,num2str(signodes(isig))); 
    end
    
end

bgdcolor = get(gcf,'color');
if ~exist('zm','var'),
    funlabel=uicontrol('style','text','units','norm','position',[.7 .95 .3 .05],...
        'string',label,'fontsize',14,'BackgroundColor',bgdcolor);
else
    if zm==0
        funlabel=uicontrol('style','text','units','norm','position',[.7 .95 .3 .05],...
            'string',label,'fontsize',14,'BackgroundColor',bgdcolor);
    else
        funlabel=uicontrol('style','text','units','norm','position',[.7 .95 .3 .05],...
            'string',['zm;' label],'fontsize',14,'BackgroundColor',bgdcolor);
    end
end


xl=[0.01 0.01  0.857 ];
yl=[1 0 0 ];
set(line(xl,yl),'color','k')

text(0.05, 0.6,sprintf('%0.3f mV',1/funscal))

if ~exist('sint') || isempty(sint), sint=1; end

text(0.9-LAY(1,1)/10, -.25 ,sprintf('%0.0f %s',nt*sint,' ms'))

if listextr~=0,
    laag=min(min(PHI));
    hoog=max(max(PHI));
    text(grid(1)/2, -.7,sprintf('%s %0.2f  to  %0.2f mV',' range ',laag,hoog))
end
axis tight
