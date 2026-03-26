% baselsplines.m
% 20121029; A. van Oosterom
% function PSI=baselsplines(PHI,zerotims,mode,fig)
% baselinecor correction using spline interpolation
% based on values at timepoints: zerotims
% if mode==1 (default) : natural spline
%    else: zero gradients enforced at the two extremes of zerotims
% splines may be monitored in figure:fig

function PSI=baselsplines(PHI,zerotims,mode,fig)

if nargin<3, mode=2; end
nt=size(PHI,2);
nl=size(PHI,1);
xx=1:nt;

if mode==1,
    % basic, natural spline: zero second derivatives at border pivots
    BASEL=spline(zerotims ,PHI(:,zerotims),xx);
else,
    % zero gradient at the border pivots
    BASEL=zeros(nl,nt);
    for i=1:nl,
        cs=spline(zerotims,[0 PHI(i,zerotims) 0]);
        BASEL(i,:)=ppval(cs,xx);
    end
end

nz=size(zerotims,2);

if nargin>3,
    
    figure(fig)
    clf
    plot(BASEL')
    hold on
    for i=1:nz,
        plot([zerotims(i) zerotims(i)],[-0.1 0.1],'k')
    end
    
    title('subtracted baselines')
end

PSI=PHI-BASEL;




% standard MATLB example
%  x = -4:4; y = [0 .15 1.12 2.36 2.36 1.46 .49 .06 0];
%        cs = spline(x,[0 y 0]);
%        xx = linspace(-4,4,101);
%        plot(x,y,'o',xx,ppval(cs,xx),'-');




