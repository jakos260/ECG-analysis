function meanVector = compareVectorLoops(varargin)
warning off;
cols=['b';'r';'k';'g';'m';'y'];

VECTORS =[];
nmap=0;
dovcg = 1;
myLegend =[];
mycolor=[];
if length(varargin) < 1
    error('This routine needs at least two parameters');
else
    pp=1;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'vcg'
                    dovcg=varargin{pp+1};pp=pp+2;
                case 'color'
                    mycolor = varargin{pp+1};pp=pp+2;
                case 'legend'
                    myLegend = varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        else
            if iscell(varargin{pp})
                
                for i=1:length(varargin{pp})
                    nmap = nmap+1; 
                    V = varargin{i};
                    V = bsxfun(@minus,V,V(1,:));
                    VECTORS{nmap} = V;
                    VECTORS{nmap} = VECTORS{nmap} / max(norm3d( VECTORS{nmap} ));                    
                end
            else
                nmap = nmap+1;               
                V = varargin{pp};
                V = bsxfun(@minus,V,V(1,:));
                VECTORS{nmap} = V;
                VECTORS{nmap} = VECTORS{nmap} / max(norm3d( VECTORS{nmap} ));               
            end           
            pp=pp+1;
        end
    end
end

meanVector = [];
for i=1:length(VECTORS)
    VECTORS{i} = bsxfun(@minus, VECTORS{i}, VECTORS{i}(1,:));
    meanVector(i,:) = sum(VECTORS{i});
    meanVector(i,:) = meanVector(i,:) / norm3d( meanVector(i,:) );
end


clf
% the horizontal plane points forward. therefore the x values have to be
% inverted for vizualization only. Then teh match the sagital plane seen
% from the head

subplot(1,3,1);
hold on;
axis off square;
title('horizontal plane');
for i=1:length(VECTORS)
    color = cols(i,:);
    if ~isempty(mycolor)
        color = mycolor(i,:);
    end
    plot(VECTORS{i}(:,2), VECTORS{i}(:,1),color,'linewidth',1);
%     plot(VECTORS{i}(1,2), -VECTORS{i}(1,1),'ko','linewidth',1,'Markersize',10);
    arrow2D([0,0],[meanVector(i,2)  meanVector(i,1)] ,1,color);
end
maxL =1;
plot([-maxL maxL],[0 0],'k');plot([0 0],[-maxL maxL],'k');axis([-maxL maxL -maxL maxL]); axis off square;

%% sagital plane

subplot(1,3,2);hold on;

title('right sagital plane');

for i=1:length(VECTORS)
    color = cols(i,:);
    if ~isempty(mycolor)
        color = mycolor(i,:);
    end
    plot(VECTORS{i}(:,1), VECTORS{i}(:,3),color,'linewidth',1);
%     plot(VECTORS{i}(1,1), VECTORS{i}(1,3),'ko','linewidth',1,'Markersize',10);
    arrow2D([0,0],[meanVector(i,1)  meanVector(i,3)] ,1,color);
end
maxL =1;
plot([-maxL maxL],[0 0],'k');plot([0 0],[-maxL maxL],'k');axis([-maxL maxL -maxL maxL]); axis off square;

%% frontal plane

subplot(1,3,3);hold on;
title('frontal plane');

for i=1:length(VECTORS)
    color = cols(i,:);
    if ~isempty(mycolor)
        color = mycolor(i,:);
    end
    plot(VECTORS{i}(:,2), VECTORS{i}(:,3),color,'linewidth',1);
end
for i=1:length(VECTORS)
    color = cols(i,:);
    if ~isempty(mycolor)
        color = mycolor(i,:);
    end
%     plot(VECTORS{i}(1,2), -VECTORS{i}(1,3),'ko','linewidth',1,'Markersize',10);
    arrow2D([0,0],[meanVector(i,2)  meanVector(i,3)] ,1,color);
end


maxL =1;
plot([-maxL maxL],[0 0],'k');plot([0 0],[-maxL maxL],'k');axis([-maxL maxL -maxL maxL]); axis off square;
if ~isempty(myLegend)
    legend(myLegend,'location','northwest','box','off');
end

%%
% plot_arrow.m
% function plot_arrow(tail,  head, linw, fracs,fillcolor);
% plot a single arrow from vectors: tail to head;  in figure: fig; linewidth: linw
% size of figure head are fractions of the length of the arrow
%      default: fracs=[0.9 0.025 0.93]; fracs(2) is adapted to axes_length ratio 
% fillcolor sets the color of the arrow head; default: 'b'
% plots in 2D supported; for 3D: see: make_arrow.m
% arrow should be contained in axis of calling script


function arrow2D(tail, head, linw, fillcolor)

fracs=[0.9 0.025 0.93];
if nargin<4,
      fillcolor='none';
end
[ni,nj]=size(tail);
if ni>nj,
   tail=tail'; 
   head=head'; 
end

vec=head-tail;
assen=axis;
lx=norm(assen(1)-assen(2));
ly=norm(assen(3)-assen(4));

xpl=vec(1)/lx; ypl=vec(2)/ly;
long=sqrt(xpl^2+ypl^2);
alpha=atan2(ypl,xpl);
normal=long*[-sin(alpha) cos(alpha)];

flank=fracs(2)*[normal(1)*lx normal(2)*ly];


ARROW=[tail; 
       tail*(1-fracs(3))+head*fracs(3);
       tail*(1-fracs(1))+head*fracs(1) + flank;
       head;
       tail*(1-fracs(1))+head*fracs(1) - flank;
       tail*(1-fracs(3))+head*fracs(3);];
   
POLYGON=ARROW(2:6,:);

if ~strcmp(fillcolor,'none'),
   fill(POLYGON(:,1),POLYGON(:,2),fillcolor);
   plot(ARROW(:,1),ARROW(:,2),'color',fillcolor,'linewidth',linw);

else
    plot(ARROW(:,1),ARROW(:,2),'color','k','linewidth',linw);

end
      
   





