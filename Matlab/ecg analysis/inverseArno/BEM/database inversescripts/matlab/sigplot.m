% sigplot voor matlab
% function extr = sigplot(PHI,label,LAY,funscal,col,leadnos,markt,listextr,linew)
%
% 20130912
%
% if leadnos==0; no signal labels plotted

function extr = sigplot(PHI,label,LAY,funscal,col,leadnos,markt,listextr,linew)
global sint;

grid        = LAY(1,:);
LAY(1,:)    = [];
nplts       = size(LAY,1);

[nsigs nt]  = size(PHI);
nplts       = min(nplts,nsigs); % note: nplts may be smaller than nlds

if nargin < 9,          linew   = 1;        end
if ~exist('label'),     label   =' ';       end
if ~exist('funscal'),   funscal = 1;        end
if ~exist('col'),       col     = 'blue';   end
if ~exist('leadnos'),   leadnos = 0;        end
if ~exist('markt'),     markt   = 0;        end
if ~exist('listextr'),  listextr= 0;        end

if length(leadnos==1),
    if leadnos ~= 0,
        if nplts == 12,
            % standard leads
            % test for augmented extrmities
            augm = 0;
            norm((PHI(5,:)-PHI(4,:)-PHI(1,:))/nt);
            
            if norm((PHI(5,:)-PHI(4,:)-PHI(1,:))/nt)>1.e-4,
                augm    = 1;
                LAB     = ['  I'; ' II'; 'III'; 'AVR'; 'AVL'; 'AVF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
            else,
                LAB     = ['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';];
            end
        end
        
        if size(LAY,1) == 9, LAB = [' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';]; end
        
        % common reference leads; OACG
        if nplts == 13, LAB = [' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';' A1';' A2';' A3';' A4'; ]; end
        
        % standard 12 signals plus four OACG leads
        if nplts == 16, LAB = ['  I'; ' II'; 'III'; ' VR'; ' VL'; ' VF'; ' V1'; ' V2'; ' V3'; ' V4'; ' V5'; ' V6';' A1';' A2';' A3';' A4'; ]; end
    end
end

% create figure with grid:
axis([0 grid(1) 0 2*grid(2)] );
axis manual
axis('off')
hold on
t           = 1:nt;
ngridcols   = grid(1);
t           = t*.8/nt*ngridcols/(ngridcols-0.2);

for i = 1:nplts, yshift(i) = 2*grid(2)+1-2*LAY(i,2); xshift(i) = LAY(i,1)-1; end

for i = 1:nplts,
    isig    = LAY(i,3);
    j       = i;
    
    plot(t+xshift(j),funscal*PHI(isig,:)+yshift(j),col,'linewidth',linew)
    plot([t(1)+xshift(j) t(nt)+xshift(j)],[yshift(j) yshift(j)],':k')
    
    if leadnos == 1,
        if nplts < 20 && nplts > 1,
            if ismember(size(LAB,1),[9 12 13 16]), text(xshift(j)+.7,yshift(j)+.25,LAB(j,:)); end
        else
            text(xshift(j)+.7,yshift(j)+.25,num2str(j));
        end
    end
    
    if length(leadnos)>1, text(xshift(j)+.7,yshift(j)+.25,num2str(isig)); end
    
    if listextr ~= 0,
        extr = extremes(PHI);
        listextr
        sprintf('Extremes: min %5.3f  at: %5d %5d; max %5.3f  at: %5d  %5d', extr)
        plot(t(extr(3))+xshift(extr(2)),funscal*extr(1)+yshift(extr(2)),'*b')
        plot(t(extr(6))+xshift(extr(5)),funscal*extr(4)+yshift(extr(5)),'*r')
    end
    
    if markt > 0 && markt <= nt,
        tmark = t(markt);
        MARKS = [];
        for ii = 1:nplts,
            j       = LAY(ii,3);
            ymark   = funscal*PHI(j,markt)+yshift(j);
            YMARK   = [ymark-.1; ymark+.1];
            TMARK   = tmark;
            TMARK   = ones(2,1)*TMARK;
            MARKS   = [MARKS plot(TMARK+xshift(j),YMARK,'r')];
        end
    end
end

funlabel = uicontrol('style','text','units','norm','position',[.7 .95 .3 .05],'string',label,'fontsize',8);

xl = [0.01 0.01  0.8];
yl = [1 0 0 ];
set(line(xl,yl),'color','k')
text(0.05, 0.6,sprintf('%0.3f mV',1/funscal))

if ~exist('sint'), sint = 1; end

text(0.5, -.25 ,sprintf('%0.0f %s',(nt-1)*sint,' [ms]'))

if listextr ~= 0,
    laag = min(min(PHI));
    hoog = max(max(PHI));
    text(grid(1)/2, -.7,sprintf('%s %0.2f  to  %0.2f mV',' range ',laag,hoog))
end

