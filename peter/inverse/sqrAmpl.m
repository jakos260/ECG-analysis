function [amplitudePotfinal,PHIbsm]=sqrAmpl(GEOM,usetimes)
% inverse determination of timing of Equivalent Double Layer source EDL; 
% optimization based on basic Marquardt type of solving thenon-linear 
% estimation problem alternate iterations between depolarization and 
% repolarization
% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

INV.lpass=5;
% versie=1;	
% PARAMETERS
maxiter=100;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;

% INV.BSM=GEOM.BSM(:,GEOM.specs(2):end);
INV.usetimes= usetimes;
INV.BSM = GEOM.BSM;
INV.ADJ=GEOM.ADJ;
if size(GEOM.BSMall,1) == size(GEOM.TVER,1)
    PHIbsm=mean(GEOM.BSMall(:,usetimes),2);
    figure(1);clf;
    showpatch(GEOM.TVER,GEOM.TITRI,PHIbsm,'range',[-0.3 0.3],loadmat('avopot.mcm'),'sublines',0.05)
    drawnow
end
INV.PHI=INV.BSM(:,usetimes);
% INV.PHI= INV.PHI(leads,:);

INV.normphi=norm(INV.PHI,'fro');
INV.T=ones(size(GEOM.AMA,2),1)*(INV.usetimes);
% compute regularization operator, e.g.surface Laplacian
[INV.REGOP]=calcREGOP(GEOM);
%precompute large matrixes
INV.AMA=GEOM.AMA;

INV.ATA=INV.AMA'*INV.AMA;  
INV.REG=INV.REGOP'*INV.REGOP;




% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=5e-6;

%% read initial estimates; 

% regularization parameters; must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error
lambopt=0.1;  
% depolarization

% Jpoint
OPT.AMPL.pot=zeros(size(GEOM.VER,1),1);
% OPT.AMPL.pot=GEOM.amplitudepot';	
OPT.AMPL.pot=0.99*ones(size(GEOM.VER,1),1);

OPT.AMPL.lambopt=lambopt;
OPT.AMPL.mu=1e-5;
OPT.AMPL.type='amplitude';

%%    
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate

%%	
% iterative approach to solving non-linear parameter estimation
% alternate between dep and rep

loopstart=1;iter=0;score=1; % start outer loop
while iter<maxiter && score
	iter=iter+1;
	[score,OPT.AMPL,TST]=optimizeNotch(INV,OPT.AMPL);				
end 
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TST=gettres_pot(INV,OPT.AMPL);
amplitudePotfinal=OPT.AMPL.pot;
figure(2);showAtria(GEOM.VER,GEOM.ITRI,amplitudePotfinal*100,'range',[70 100])
rdfinal=TST.rd;
COR=corrcoef(TST.PHIA,INV.PHI);corfinal=COR(2,1);
disp(num2str([ iter rdfinal corfinal]))

function [score,amplitudeopt,TST]=optimizeNotch(INV,amplitudeopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,amplitudeopt);
starttres=TST.tres;
score=0;
lamb=amplitudeopt.lambopt; 

Splus=getS(INV.T,amplitudeopt.pot+0.01);
Smin =getS(INV.T,amplitudeopt.pot-0.01);
Sprime=(Splus-Smin)/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+amplitudeopt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-amplitudeopt.mu^2*INV.REG*amplitudeopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	delamplitude=inv(GTGML)*righths;
% 	weight function to limit plateau phase apmlitue to at least 85%
	newpot=amplitudeopt.pot+delamplitude; 
	newpot(newpot<0)=0;
	newpot(newpot>1)=1;
	testopt=amplitudeopt;testopt.pot=newpot;
	TST=gettres_pot(INV,testopt);
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			amplitudeopt.pot=newpot; 
			score=1;
		end
		break;
	end
end  % end of inner loop
amplitudeopt.lambopt=lamb/4;
%% ======================================================================
function TST=gettres_pot(INV,opt)


S=getS(INV.T,opt.pot);
TST.PHIA = INV.AMA*S;
TST.RES = INV.PHI - TST.PHIA;
TST.rd=norm(TST.RES,'fro')/norm(INV.PHI,'fro'); %NOTE: unfiltered rd
TST.reg=norm(INV.REGOP*opt.pot)/1000;
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;
%% =======================================================================
function REGOP=calcREGOP(GEOM)

if max(max(abs(GEOM.VER)))>1
	REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
else
	REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
end

%% =======================================================================
   
function S=getS(T,amplitudepot)

S=ones(size(T));
for i=1:size(S,1)
	S(i,:)=S(i,:)*amplitudepot(i);
end


