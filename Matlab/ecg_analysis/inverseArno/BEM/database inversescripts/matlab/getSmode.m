function S = getSmode(T,dep,rep,SPECS,mode,varargin)
%
% INPUT: 
% T     = matrix of # vertices * timescale QRT
% dep   = depolarization times
% rep   = repolarization times
% SPECS = 
% mode  = use only depolarization or also repolarization
%
% variant of gets dedicated to ventricles
% function S = getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11
%
% pS(3);    % INPUT determines the slope leading up to the apex
% pS(4);    % INPUT determines the (negative) slope following the apex

nt = size(T,2);
if size(dep,1) < size(dep,2), dep = dep'; end
if size(rep,1) < size(rep,2), rep = rep'; end

%% depolarization

% Anti-Hakkel: voorlopig niet gebruiken!!
% if length(varargin) == 1, % if 6 inputs are given instead of 5.
%     GEOM    = varargin{1};
%     S       = getSpc(T,dep,rep,SPECS,mode,GEOM.neigh);
%     % S       = getSmodeHakkel(T,dep,rep,SPECS,mode,varargin{1});
%     10
%     return; 
% end

TDEP = T - dep*ones(1,nt);
if size(SPECS.depSlope,1) == 1,
    S = 1./(1+exp(-SPECS.depSlope.*TDEP));
else
    S = 1./(1+exp(bsxfun(@times,TDEP,-SPECS.depSlope(:))));
end

%% repolarization
if mode >= 4,
    TREP = T-rep*ones(1,nt)-SPECS.repCorrection;
    if SPECS.useCumsum,
        Y = (SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP)); 
        Y = 1 - bsxfun(@times,cumsum(Y,2),1./sum(Y,2));
    else
        Y = -(SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP)); 
        Y = bsxfun(@times,Y,1./Y(:,1));
    end
    
    % apply the up-slope: logistic shape;   
    S = S .* Y;
    if ~isempty(SPECS.scaleAmpl),
        S = bsxfun(@rdivide,S,max(S,[],2)./SPECS.scaleAmpl);   % oostep1: 20121015 much faster
    else
        S = bsxfun(@rdivide,S,max(S,[],2));                    % oostep1: 20121015 much faster
    end
end