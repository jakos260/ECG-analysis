% plot_arrow.m
% function plot_arrow(tail,  head, fig , linw, fracs,fillcolor);
% plot a single arrow from vectors: tail to head;  in figure: fig; linewidth: linw
% size of figure head are fractions of the length of the arrow
%      default: fracs=[0.9 0.025 0.93]; fracs(2) is adapted to axes_length ratio 
% fillcolor sets the color of the arrow head; default: 'b'
% plots in 2D supported; for 3D: see: make_arrow.m
% arrow should be contained in axis of calling script

% 2011-06-01; A. van Oosterom

function plot_arrow(tail, head, fig, linw, fracs, fillcolor);

if nargin<5,
   fracs=[0.9 0.025 0.93];
end
if isempty(fracs),
   fracs=[0.9 0.025 0.93];
end

figure(fig)
hold on
[ni,nj]=size(tail);
if ni>nj,
   tail=tail'; head='head'; 
   [ni,nj]=size(tail);
end
dim=nj;

vec=head-tail;
long=norm(vec);
assen=axis
lx=norm(assen(1)-assen(2))
ly=norm(assen(3)-assen(4))

xpl=vec(1)/lx; ypl=vec(2)/ly;
long=sqrt(xpl^2+ypl^2);
alpha=atan2(ypl,xpl);
normal=long*[-sin(alpha) cos(alpha)];

flank=fracs(2)*[normal(1)*lx normal(2)*ly];

ARROW=[tail; 
       head;
       tail*(1-fracs(1))+head*fracs(1) + flank;
       tail*(1-fracs(3))+head*fracs(3);
       tail*(1-fracs(1))+head*fracs(1) - flank;
       head;];
POLYGON=ARROW(2:6,:);
plot(ARROW(:,1),ARROW(:,2),'linewidth',linw);
if nargin<6,
      fillcolor='b';
end
if ~strcmp(fillcolor,'none'),
   fill(POLYGON(:,1),POLYGON(:,2),fillcolor);
end
      
   





