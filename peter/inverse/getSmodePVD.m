function S=getSmode(T,dep,rep,pS,scaleAmpl,mode,varargin)
% variant of gets dedicated to ventricles
% function S=getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11

% global GEOM



pup= -2;	% determines the slope of the upstroke

nt=size(T,2);
[ni nj]=size(dep);if ni<nj, dep=dep';end
[ni nj]=size(rep);if ni<nj, rep=rep';end

%% depolarization

TDEP=T-dep*ones(1,nt);
if max(size(pup))==1
    S=1./(1+exp(pup.*TDEP));
else
    S=1./(1+exp(bsxfun(@times,TDEP,pup)));
end

if mode==1
    % 	S=diag(1./max(S,[],2))*S;
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
end
%% repolarization
% down
if mode>=4
    pdp(2)=pS(1);    % INPUT determines the slope leading up to the apex
    pdp(3)=pS(2);    % INPUT determines the (negative) slope following the apex
    TREP=T-rep*ones(1,nt);
    Y= (1./(1+exp(pdp(2)*TREP))).*1./(1+exp(pdp(3)*TREP));
    
    if mode==6
        D=bsxfun(@times,TREP,(-0.3.*max(Y,[],2)) ./ (rep-dep)); % oostep1: 20121015 much faster
        D(D<0)=0;
        Y=D+Y;
    end
    %     Y= (1./(1+exp(p3*TREP)))./(1+exp(p4*TREP));
    
    S=S.*Y;
    if ~isempty(scaleAmpl)
        %         S=diag(scaleAmpl./max(S,[],2))*S;
        S=bsxfun(@rdivide,S,max(S,[],2)./scaleAmpl); % oostep1: 20121015 much faster
    else
        %         S=diag(1./max(S,[],2))*S;
        S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
    end
end

