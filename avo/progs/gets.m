% gets.m
% variant of gets
% function S=gets(T,dep,rep,p,mode)
% t: rowvector; dep and rep: columnvectors
% A. van Oosterom; 2017-03-08;
%
% bsxfun(@rdivide, .. introduced at line 61 
% for fine tuning p(5) assign a negative sign to p(5) in calling script


function S=gets(T,dep,rep,p,mode)
nt=size(T,2);
nn=size(T,1);

if size(dep,1)~=nn, dep=dep'; end
if size(rep,1)~=nn, rep=rep'; end

TDEP=T-dep*ones(1,nt);

% p(1): slope upstroke*4;
% p(2): init value Y_dom, setting initial downslope of TMP
% p(3): parameter setting leading slope of Tdom;
% p(4): parameter setting trailing slope of Tdom;
% p(5): extra shift forcing the timing of the apex to coincide toward that of
%       apex Tdom
% if size(p,1)>5
% p(6): amplitude U-wave (Gauss)
% p(7): width (sd) of U wave
% p(8): timing of peak U_wave relative to rep

if mode~=1,

    TREP=T-rep*ones(1,nt)-p(5);
    % compute (-1*) derivative of the downward slope of the TMP;
    
    Y=(p(2)+1./(1+exp(p(3)*TREP)))./(1+exp(p(4)*TREP));
    
    % if size(p,1)>5,  
    %     TU=T-p(8)*ones(nn,nt);   
    %     Uwave=p(6)*exp(-0.5*(TU/p(7)).^2);
    %     Y=Y+Uwave;
    % end
    
    % compute TMP; unit below the curve
    % Y=(1-diag(1./sum(Y,2))*cumsum(Y')');
    
    % CS=cumsum(Y')';
    % Y=1-bsxfun(@rdivide,CS,max(CS,[],2));
    
    % apply the up-slope: logistic shape;
    S=(1-Y)./(1+exp(-p(1)*TDEP));
   
    % force a unit upstroke
    S=bsxfun(@rdivide,S,max(S,[],2));   
else
    S=1./(1+exp(-p(1)*TDEP));
end








% prullies

% (re-)establish unit upstroke
% S=diag(1./max(S'))*S;
 
%     if sign(p(5))<0,
%         p(5)=abs(p(5));
%         % recompute p(5)
%         mean_rep=mean(rep);
%         trep=T(1,:)-(mean_rep+p(5))*ones(1,nt);
%         y=(p(2)+1./(1+exp(p(3)*trep)))./(1+exp(p(4)*trep));
%         [ma,ima]=max(y);
%         p(5)=p(5)-ima+mean_rep;
%     end
% pause
    
