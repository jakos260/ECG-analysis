% plot12lds.m
% function plot12(PHI,label,mode,funscal,color,linw,fonts);
% assumes inclusion in PHI:  I,II,III,VRA,VLA,VLF,V1, ...V6
% if mode=1, minimap format

function plot12lds(PHI,label,mode,funscal,color,linw,fonts);
if nargin<7, fonts=10; end
if nargin<6, linw=1; end
if nargin<5, color='b'; end
if nargin<4, funscal=1; end
if nargin<3, mode=0; end

IDENT=['  I';' II';'III';'VRA';'VLA';'VLF';' V1';' V2';' V3';' V4';' V5';' V6'];

% prepare plotting 12 leads
[nl ntpl]=size(PHI);
tt=1:ntpl;
tt=tt*.8/ntpl;

% ax2= axes('units','norm','pos',[0.1 0.2 0.8 0.7]);
% axes(ax2)
axis([0 8 0 6])

axis manual
axis off
hold on

title(label,'fontsize', fonts)

if mode==1,
   % set the data required for plotting the 12 leads in minimap format
   yshift=[ 5.5 1    3.5 5.5 5.5 0.2   4  4  3.25 2.5 2.5 2.5 ];
   xshift=[ 3.5 3.5  7   0   7   7     1  2  3    4   5   6  ];
end

if mode~=1,
   % set the data required for plotting the 12 leads in standard format
   yshift=[5 3 1  5 3 1  5 3 1  5 3 1];
   xshift=[1 1 1  3 3 3  5 5 5  7 7 7];
end


XSHIFT=xshift'*ones(1,ntpl);
YSHIFT=yshift'*ones(1,ntpl);
T=ones(12,1)*tt;
TDAT12=XSHIFT+T;
basel=plot(([XSHIFT, XSHIFT+0.8])',([YSHIFT YSHIFT])','k:');
txtime=text(0.45,0.8,[num2str(ntpl) 'ms']);
ampli=plot([0 0 0.85],[2 1 1],'k'); 
txtmv=text(0.2,2,[num2str(1/funscal) 'mV']);
sigs=(funscal*PHI+YSHIFT)';
plts=plot((TDAT12)',sigs,color, 'linewidth',linw);
if fonts>0, txtlds=text(xshift+0.7,yshift +.5,IDENT(:,1:3),'fontsize',fonts);, end;
