function S=getSvplateau(T,dep,rep,pS,scaleAmpl,mode)
% variant of gets dedicated to ventricles
% function S=getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11
global dodeprep
pup(1)= -2;	% determines the slope of the upstroke 


%plateau
plp(1)=0.25;
plp(2)=100;
nt=size(T,2);
[ni nj]=size(dep);if ni<nj, dep=dep';end
[ni nj]=size(rep);if ni<nj, rep=rep';end


%% depolarization
TDEP=T-dep*ones(1,nt);
S=1./(1+exp(pup(1)*TDEP));	
if mode==1	
	amplitude=max(S,[],2);
	S=diag(1./amplitude)*S;
	return;
end
%% plateau
if exist('dodeprep','var') && dodeprep==0
	plp(1)=.08;
	plp(2)=8;
	TPLATEAU=T-(dep+plp(2))*ones(1,nt);	
	P=1-plp(1)./(1+exp(-.35*TPLATEAU)); 
% 	N= (1./(1+exp(-0.15*(TDEP-25))))+(1.0./(1+exp(0.025*(TDEP-60))));
% 	N= N-1;
% 	amplitude=max(N,[],2); 
% 	N=1+0.05*diag(1./amplitude)*N;
	S=S.*P;
end

%% repolarization
% down
pdp(2)=pS(1);    % INPUT determines the slope leading up to the apex
pdp(3)=pS(2);    % INPUT determines the (negative) slope following the apex

TREP=T-rep*ones(1,nt); 
if exist('dodeprep','var') && dodeprep==1
	Y= (1./(1+exp(pdp(2)*TREP))).*1./(1+exp(pdp(3)*TREP));
else
	Y= (1./(1+exp(.5*pdp(2)*TREP))).*1./(1+exp(1*pdp(3)*TREP));
end
S=S.*Y;
amplitude=max(S,[],2); 
S=diag(1./amplitude)*S;
% if size(T,1)==1, plot(Y,'g','linewidth',2);hold on; plot(P,'r','linewidth',2);end


if ~isempty(scaleAmpl)
	if max(scaleAmpl)>1
		S=diag(scaleAmpl/100)*S;
		S=S+(1-scaleAmpl/100)*ones(1,size(S,2));
	else
		S=diag(scaleAmpl)*S;
		S=S+(1-scaleAmpl)*ones(1,size(S,2));
	end

end
