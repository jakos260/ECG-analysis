% splines_plus.m
% 20121028; A. van Oosterom
% function F_int=splines_plus(F,pivots,xx,mode)
% if mode==1 (default) : natural spline
%    else: zero gradients enforced at the two extremes of pivots

function F_int=splines_plus(F,pivots,xx,mode)

% 
% function PSI=baselsplines(PHI,zerotims,mode)
% 
% if nargin<3, mode=2; end
% nt=size(PHI,2);
% nl=size(PHI,1);
% xx=1:nt;
% 
% if mode==1,
%     % basic, natural spline
%     BASEL=spline(zerotims ,PHI(:,zerotims),xx);
% else,
%     % spline using zero gradient at the bordering pivots 
%     BASEL=zeros(nl,nt);
%     for i=1:nl,
%         cs=spline(zerotims,[0 PHI(i,zerotims) 0]);
%         BASEL(i,:)=ppval(cs,xx);
%     end
% end
% 


% F=F';
% pivots
% xx

if nargin<3, mode=1; end
nj=size(F,2);
ni=size(F,1);
nn=max(size(xx));


if mode==1,
    % basic, natural spline
    F_int=spline(pivots,F,xx);
else,
    % spline including minimization of gradient at the border pivots 
    F_int=zeros(ni,nn);
    for i=1:ni,
        cs=spline(pivots,[0 F(i,:)  0]);
        F_int(i,:)=ppval(cs,xx);
    end
end
