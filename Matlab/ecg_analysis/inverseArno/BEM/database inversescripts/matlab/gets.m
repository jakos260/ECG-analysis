% gets.m
% variant of gets 
% function S=gets(t,dep,rep,p,mode)
% t: rowvector; dep and rep: columnvectors
% A. van Oosterom; 2008-12-04; 
% 2012-10_23; speeded up by PO

function S = gets(T,dep,rep,p,mode)

nt = size(T,2);
nn = size(T,1);

if size(dep,1)<size(dep,2), dep = dep'; end
if size(rep,1)<size(rep,2), rep = rep'; end


TDEP = T-dep*ones(1,nt);

% p(1): slope upstroke*4;
% p(2): init value Y_dom, setting initial downslope of TMP
% p(3): parameter setting leading curvature of Tdom;
% p(4): parameter setting trailing curvature of Tdom;
% p(5): extra shift forcing the timing of the apex to coincide toward that of 
%       apex Tdom 
   
if mode ~= 1,
    TREP    = T-rep*ones(1,nt)-p(5); 
    p3      = ones(nn,1)*p(3); 
    p4      = ones(nn,1)*p(4);

    % compute (-1*) derivative of the downward slope of the TMP; 
    Y       = (p(2)+1./(1+exp(diag(p3)*TREP)))./(1+exp(diag(p4)*TREP)); 

    % compute TMP; unit upstroke
    Y       = (1-diag(1./sum(Y'))*cumsum(Y')');
    % apply the up-slope: logistic shape;  

    S       = Y./(1+exp(-p(1)*TDEP));
    % re-establish unit upstroke
    % S = diag(1./max(S'))*S;
    S       = bsxfun(@rdivide,S,max(S,[],2));
    
else
    S       = 1./(1+exp(-p(1)*TDEP));
end

