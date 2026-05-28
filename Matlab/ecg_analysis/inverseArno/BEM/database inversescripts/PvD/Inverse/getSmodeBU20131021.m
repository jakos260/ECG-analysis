function S=getSmode(T,dep,rep,pS,scaleAmpl,mode)
% variant of gets dedicated to ventricles
% function S=getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11
pup(1)= -4;	% determines the slope of the upstroke 
%plateau

nt=size(T,2);
[ni nj]=size(dep);if ni<nj, dep=dep';end
[ni nj]=size(rep);if ni<nj, rep=rep';end


%% depolarization
% TDEP=T-dep*ones(1,nt);
TDEP=bsxfun(@minus,T,dep); % oostep1 20121016
S=1./(1+exp(pup(1)*TDEP));




if mode==1	
% 	S=diag(1./max(S,[],2))*S;
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster  
end
%% repolarization
% down
if mode>=4
    pdp(2)=pS(1);    % INPUT determines the slope leading up to the apex
    pdp(3)=pS(2);    % INPUT determines the (negative) slope following the apex
%     TREP=T-rep*ones(1,nt); 
    TREP=bsxfun(@minus,T,rep); % oostep1 20121016 faster

    Y= (1./(1+exp(pdp(2)*TREP))).*1./(1+exp(pdp(3)*TREP));
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

return
%% plateau
if mode == 4
    if ~isempty(scaleAmpl)
        S=S.*(scaleAmpl * ones(1,size(S,2)));
    end
elseif mode==2
	plp=scaleAmpl;
	TPLATEAU=T-(dep+6)*ones(1,nt);	
	P=(plp*ones(1,nt))./(1+exp(-1*TPLATEAU)); 
	temp=max(S.*(1-P),[],2);
	for i=1:size(P,1)
		P(i,:)=1-(P(i,:)/temp(i));
	end
	S=S.*P;
	
% 	S=diag(1./max(S,[],2))*S;
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster

    
	return
elseif mode==3
	N= (1./(1+exp(-0.05*(TDEP-25))))+(1.0./(1+exp(0.025*(TDEP-60))));
	N= N-1;
	N=1+scaleAmpl*diag(1./max(N,[],2))*N;
	S=S.*N;
% 	S=diag(1./max(S,[],2))*S;	
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster

% 	return
elseif mode==5
	TPLATEAU=T-(dep+5)*ones(1,nt);	
	P=1-(scaleAmpl*ones(1,nt))./(1+exp(-1*TPLATEAU)); 
	S=S.*P;
% 	S=diag(1./max(S,[],2))*S;
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster

	
	scaleAmpl=max(scaleAmpl*100/(rep-dep),0);
	N= (1./(1+exp(-0.15*(TDEP-25))))+(1.0./(1+exp(0.03*(TDEP-50))));
	N= N-1;
	N=1+scaleAmpl*diag(1./max(N,[],2))*N;
 	S=S.*N;
% 	S=diag(1./max(S,[],2))*S;	
    S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster

end

	

	

% if size(T,1)==1, plot(Y,'g','linewidth',2);hold on; plot(P,'r','linewidth',2);end

% if exist('dodeprep','var') && dodeprep==1 && ~isempty(scaleAmpl)
% 	for i=1:size(S,1)
% 		S(i,:)=S(i,:)*scaleAmpl(i);
% 	end
% end

% S=S-0.85;