function [depfinal,repfinal,notchPotfinal,corfinal,rdfinal,iterfinal,RESNOW]=ecginv_atria(GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
% inverse determination of timing of Equivalent Double Layer source EDL; 
% optimization based on basic Marquardt type of solving thenon-linear 
% estimation problem alternate iterations between depolarization and 
% repolarization
% peter van Dam; 20080929

% if INV.mode==1, loopstart=1;nloop=1; repscore=0;end % dep only
% if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
% if INV.mode==3, loopstart=3;nloop=3; end % notch
% if INV.mode==4, loopstart=1;nloop=2; end % dep rep 
% if INV.mode==5, loopstart=1;nloop=3; end % dep rep & notch
% if INV.mode==6, loopstart=6;nloop=6; end % AP amplitude only
global NFIG
global lpass	
global dodeprep

versie=1;	
casedir	='invresults\';
domonitor=0;


% PARAMETERS
lpass=5;				% moving avarage lowpass filtering
sigscal=GEOM.specs(1);				% for SIGPLOT
useSurfLapl=1;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;

INV.mode=mode;
INV.BSM=zeromean(GEOM.BSM(:,GEOM.specs(2):end));
INV.usetimes=min(usetimes,size(INV.BSM,2));
INV.PHIREF=zeromean(INV.BSM(leads,1:INV.usetimes));

INV.ADJ=GEOM.ADJ;

INV.normphi=norm(INV.PHIREF,'fro');
INV.T=ones(size(GEOM.AMA,2),1)*([1:INV.usetimes]-1);
% compute regularization operator, e.g.surface Laplacian
[INV.REGOP,INV.REGOPREP]=calcREGOP(GEOM,useSurfLapl);
%precompute large matrixes
INV.AMA=zeromean(GEOM.AMA(leads,:));
if 0%useWCT
    wct=GEOM.wct;
    INV.BSM=INV.BSM-ones(size(INV.BSM,1),1)*mean(INV.BSM(wct,:));
    INV.PHIREF=INV.PHIREF-ones(size(INV.PHIREF,1),1)*mean(INV.PHIREF(wct,:));
    INV.AMA=INV.AMA-ones(size(INV.AMA,1),1)*mean(INV.AMA(wct,:));
end
INV.ATA=INV.AMA'*INV.AMA;  
INV.REG=INV.REGOP'*INV.REGOP;
INV.REGREP=INV.REGOPREP'*INV.REGOPREP;
INV.pSn=GEOM.pS;
% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=2e-4;
MINRD=0.16;
if strcmp(GEOM.type,'_atria')
    MINRD=0.201;
end
maxiter=15; 
%% read initial estimates; 

% regularization parameters; must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error
lambopt=0.1;  
% depolarization
OPT.DEP.mu=mudep;%1.5e-6;   
% if strcmp(GEOM.type,'_atria'),OPT.DEP.mu=1.5e-4; end
OPT.DEP.lambopt=lambopt;  
OPT.DEP.tims=inittimes(:,1);
OPT.DEP.type='dep';

% repolarization
OPT.REP.tims=inittimes(:,2);
OPT.REP.lambopt=lambopt;  
OPT.REP.type='apd';
OPT.REP.mu=murep;%3.0e-5; 

if strcmp(GEOM.type,'_atria')
%     OPT.REP.mu=1.5e-4; 
end
% Jpoint
OPT.NOT.pot=zeros(size(GEOM.VER,1),1);
if INV.mode==3 || INV.mode==5 || INV.mode==6
	OPT.NOT.pot=0.01*ones(size(GEOM.VER,1),1);
end
if dodeprep==0
	ampl=sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
	OPT.NOT.pot=1-ampl;
else
	OPT.AMP.pot=zeros(size(GEOM.VER,1),1);
end
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
OPT.NOT.usetimes=qrsduration+1:qrsduration+100;
OPT.NOT.lambopt=lambopt;
OPT.NOT.mu=1e-6;
OPT.NOT.type='notch';

% AP amplitude
if dodeprep==1
OPT.AMP.pot=0.1*ones(size(GEOM.VER,1),1);
else
	if exist('ampl')
		OPT.AMP.pot=ampl;
	else
		OPT.AMP.pot=sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
	end
end
if INV.mode==7,OPT.AMP.pot=OPT.AMP.pot-2;end
OPT.AMP.usetimes=qrsduration+1:qrsduration+100;
OPT.AMP.lambopt=lambopt;
OPT.AMP.mu=1e-6;
OPT.AMP.type='amplitude';

%%
OPTSTART=OPT.DEP;
TST=gettres_v(INV,OPTSTART,OPT.REP,OPT.NOT,OPT.AMP);
%% =======================================================================
% DOCUMENT CASE  

logfile=[casedir runcase '.log'];
outfil= [casedir runcase '.src'];
outfilP= [casedir runcase '.pS'];
fh=fopen(logfile,'w');
fprintf(fh,'%s\n',['file: ' logfile]);
fprintf(fh,'%s\n',['date: ' datestr(clock)]);
fprintf(fh,'%s\n',['prog: invgeom (' num2str(versie) ')' ]);
fprintf(fh,'%s\n',['initvals:   ' initcase]);
fprintf(fh,'%s\n',['subject:   ' GEOM.subject]);

if useSurfLapl==1
	fprintf(fh,'%s\n','SL regularization');
else 
	fprintf(fh,'%s\n','inverse 1.distance^2 regularization ');
end
fprintf(fh,'mudep: %0.4g    murep: %0.4g \n',OPT.DEP.mu,OPT.REP.mu);
fprintf(fh,'%s\n','ifix(1) tfix(1) dampd   dampr usetimes   tblc   mode  maxiter');

RESNOW=[];
%%    
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate

fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
temp=[0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));

fprintf(fh,'%s\n','iter   min   range  max    std  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
txt=sprintf('%s\n','iter   min   range  max    STD  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
disp(txt);
temp=[0 min(OPT.DEP.tims) diff(range((OPT.DEP.tims))) max(OPT.DEP.tims) std(OPT.DEP.tims)...
				 min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims)) ...
				 max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
				 TST.reg TST.rd TST.tres];
txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f Start time: %s',temp,datestr(now,'HH,MM.SS'));
disp(txt);
fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
fprintf(fh,'time run: %s\n',datestr(now,'HH,MM.SS')) ;
RESNOW(1,:)=[temp norm(TST.RES)];
ik=2;
  

%% ingredients for monitor
VER=GEOM.VER;ITRI=GEOM.ITRI;
zin=trisense(GEOM.ITRI,GEOM.VER);
RIM=getrim(GEOM.ITRI,zin,[1 size(GEOM.ITRI,2)]);
% RANGES lcavity rcavity and epicardial triangle indices
nitri=length(GEOM.ITRI);	ind=1:nitri;mask=zeros(nitri,1);
range1=findendoTRi(GEOM.VER,GEOM.ITRI,GEOM.LVER);
range2=findendoTRi(GEOM.VER,GEOM.ITRI,GEOM.RVER);
for i=1:nitri,if isempty(find(i==range1)|find(i==range2)),mask(i)=1;end;end
range3=ind(mask==1);range4=1:nitri;range5=1:nitri;range6=1:nitri;	
if domonitor==1
	figure(3);clf;sigplot_p(INV.PHIREF,GEOM.subject,GEOM.LAY,sigscal,'b',1,0,1);
	sigplot_p(TST.PHIA(:,1:INV.usetimes),[initcase '  rd: ' num2str(TST.rd)],GEOM.LAY,sigscal,'r',0,0,0);	
	figure(1);clf; montims=OPT.DEP.tims; monitor;
	subplot(2,3,1); title([initcase ' producing: ' runcase '  iter:' num2str(0) ' rd:' num2str(TST.rd)])
	figure(2);clf; montims=OPT.REP.tims;   monitor;
	subplot(2,3,1); title([initcase ' producing: ' runcase '  iter:' num2str(0) ' rd:' num2str(TST.rd)])
% 	disp('arrange windows/figures'); 	pause
end
%%	
 
% iterative approach to solving non-linear parameter estimation
% alternate between dep and rep

iter=0;score=1; % start outer loop
startTime=clock;
depscore=1;repscore=1;
if INV.mode==1, loopstart=1;nloop=1; repscore=0;end	% dep only
if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
if INV.mode==3, loopstart=3;nloop=3; end			% notch only
if INV.mode==4, loopstart=1;nloop=2; end			% dep & rep 
if INV.mode==5, loopstart=1;nloop=2; end			% dep & notch & rep 
if INV.mode==6, loopstart=1;nloop=1; repscore=0;end % dep & notch
if INV.mode==7, loopstart=7;nloop=7; end			% AP amplitude only
totscore=0;
TST.rd=1;
while iter<maxiter && (depscore || repscore) %&& TST.rd>MINRD 
	iter=iter+1;
	for loop=loopstart:nloop
		switch loop
			case 1 % DEP
				[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
				if INV.mode==5 || INV.mode==6
					if dodeprep
					[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
					else
					[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
					end
				end
				opt=OPT.DEP;		
				totscore=totscore+depscore;
			case 2 %REP
				if INV.mode==8
					[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
				end		
				[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);
				opt=OPT.REP;
				totscore=totscore+repscore;				
			case 3
				[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);
				opt=OPT.DEP;															
			case 7
				[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
				[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
				[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);				
				opt=OPT.DEP;		
				totscore=totscore+depscore;
		end
		temp=[iter min(opt.tims) diff(range(opt.tims)) max(opt.tims) std(opt.tims)...
			  min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims))...
			  max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
			  TST.reg TST.rd TST.tres];
		if score==0
			k=rem(ik,2)+1;s
			temp(9:10)=RESNOW(ik-k,9:10);
			temp(11)=RESNOW(ik-1,11);
		end
		txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f run %s time %s',...
					temp,opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
		disp(txt);
		fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
		fprintf(fh,'run %s time %s\n',opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
		RESNOW(ik,:)=[temp norm(TST.RES)];
		ik=ik+1;
		if domonitor==1
			figure(loop);clf;montims=opt.tims; monitor;
			subplot(2,3,1)
			title([initcase ' producing: ' runcase '  iter:' num2str(iter) ' rd:' num2str(TST.rd)])
			figure(3);   clf
			sigplot_p(INV.PHIREF,runcase,GEOM.LAY,sigscal,'b',1,0,1);
			sigplot_p(TST.PHIA,runcase,GEOM.LAY,sigscal,'r',0,0,0);
			title('blue: measured; red: simulated');
			uicontrol('style','text','units','norm','position',[0 .95 .3 .05],... 
			'string',['rd:' num2str(TST.rd)],'fontsize',14);
			figure(NFIG+17);subplot(2,1,2);
			cols=['.k';'xr';'db';'sg';'+m';'oc'];
            if ~strfind(GEOM.type,'atria')
                for i=1:6,plot(OPT.DEP.tims(GEOM.part==i),OPT.REP.tims(GEOM.part==i)-OPT.DEP.tims(GEOM.part==i),cols(i,:));hold on;end
                legend(GEOM.parts)
            end
		end
	end
end 
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TST=gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
depfinal=OPT.DEP.tims;
repfinal=OPT.REP.tims;
% if INV.mode<7
if dodeprep
	notchPotfinal=OPT.AMP.pot;
else
	notchPotfinal=OPT.NOT.pot;
end
rdfinal=TST.rd;
COR=corrcoef(TST.PHIA,INV.PHIREF);
corfinal=COR(2,1);
iterfinal=iter;


saveasci(outfil,[depfinal repfinal]);
saveasci(outfilP,INV.pSn);
fclose(fh);
% end of the inverse proceedure; view/and store results

if 0 && domonitor
	figure(NFIG+7);HEADLN=2;
	plot(1:length(RESNOW(:,11))-HEADLN,RESNOW(HEADLN+1:end,11));grid
	figure(NFIG+3);clf
	sigplot_p(INV.PHIREF,runcase,GEOM.LAY,sigscal,'b',1,0,1);
	sigplot_p(TST.PHIA,runcase,GEOM.LAY,sigscal,'r',0,0,0);
	title('blue: measured; red: simulated');
	uicontrol('style','text','units','norm','position',[0 .95 .3 .05],... 
	'string',['rd:' num2str(TST.rd)],'fontsize',14);
	figure(NFIG+4);clf; % butterfly plot of all recorded and simulated leads
	subplot(2,1,1);plot(INV.PHIREF');title('measured');
	subplot(2,1,2);plot(TST.PHIA');title('simulated')
	figure(NFIG+8);plot(TST.S');
	figure(NFIG+6);showAtria2(GEOM.VER,GEOM.ITRI,depfinal,repfinal,'do2',1);
	if INV.mode==1
% 		plot(OPT.DEP.tims,OPT.REP.tims-OPT.DEP.tims,'.b');
	else
		figure(NFIG+17);subplot(2,1,2);
		cols=['pk';'xr';'db';'sg';'+m';'oc'];
% 		for i=1:6,plot(OPT.DEP.tims(GEOM.part==i),OPT.REP.tims(GEOM.part==i)-OPT.DEP.tims(GEOM.part==i),cols(i,:),'markersize',4);hold on;end
% 		legend(GEOM.parts);		
		plot(OPT.DEP.tims,OPT.REP.tims-OPT.DEP.tims,'.r');
	end
	if INV.mode>2
		if INV.mode<7
			S=TST.S;
			dS=diffrows(S);	
			ndepth=zeros(size(S,1),1);
			for i=1:size(S,1)
				upi=find(dS(i,:)==max(dS(i,:)));
				ndepth(i)=1- min(S(i,upi+4:upi+50));				
% 				ndepth(i)=max(S(i,upi+4:end))- min(S(i,upi+4:upi+50));
			end			
			figure(NFIG+16);showAtria2(GEOM.VER,GEOM.ITRI,repfinal-depfinal,ndepth,'do2',1);
		else
			figure(NFIG+16);showAtria2(GEOM.VER,GEOM.ITRI,repfinal-depfinal,OPT.AMP.pot,'do2',1);
		end
		pause(5)
	end
end

%% =======================================================================
function [score,amplopt,TST]=optimizeAmpl(INV,amplopt,depopt,repopt,notchopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,amplopt,depopt,repopt,notchopt);
starttres=TST.tres;
score=0;
lamb=amplopt.lambopt; 

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot,amplopt.pot+.001,4);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot,amplopt.pot-.001,4);
Sprime=(Splus(:,amplopt.usetimes)-Smin(:,amplopt.usetimes))/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+amplopt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-amplopt.mu^2*INV.REG*amplopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres
while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	deltau=inv(GTGML)*righths;
% 	weight function to limit plateau phase apmlitue to at least 85%
	newpot=amplopt.pot+deltau;
	newpot(newpot<0)=0;
	newpot(newpot>1)=1;
	testopt=amplopt;test.pot=newpot;
	TST=gettres_pot(INV,testopt,depopt,repopt,notchopt);
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			ampl.pot=newpot; 
			score=1;
			disp('optimize amplitude')			
		end
		break;
	end
end  % end of inner loop
amplopt.lambopt=lamb/4;
TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);
%% =======================================================================
function [score,notchopt,TST]=optimizeNotch(INV,notchopt,depopt,repopt,amplopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,notchopt,depopt,repopt,amplopt);
starttres=TST.tres;
score=0;
lamb=notchopt.lambopt; 

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot+0.01,amplopt.pot,4);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot-0.01,amplopt.pot,4);
Sprime=(Splus(:,notchopt.usetimes)-Smin(:,notchopt.usetimes))/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+notchopt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-notchopt.mu^2*INV.REG*notchopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	delnotch=inv(GTGML)*righths;
% 	weight function to limit plateau phase apmlitue to at least 85%
	newpot=notchopt.pot+delnotch; 
% 10 is appoximates the (liniear?) relation between notchpotential 
% 	and the parameter regulating the notch amplitude
	
	newpot(newpot<0)=0;
	newpot(newpot>1)=1;
	testopt=notchopt;testopt.pot=newpot;
	TST=gettres_pot(INV,testopt,depopt,repopt,amplopt);
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			notchopt.pot=newpot; 
			score=1;
			disp('optimize notch')
		end
		break;
	end
end  % end of inner loop
notchopt.lambopt=lamb/4;
TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);

%% =======================================================================
function [score,opt,TST]=optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration

TST=gettres_v(INV,opt,keepopt,notchopt,amplopt);
starttres=TST.tres;
score=0;
lamb=opt.lambopt; 

if strcmp(opt.type,'dep')
	Splus=getS(INV.T,opt.tims+1,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	Smin =getS(INV.T,opt.tims-1,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	startopt=opt;	
elseif strcmp(opt.type,'rep')
	Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	startopt=opt;
else %if strcmp(opt.type,'apd')
	Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
	startopt=opt;
	startopt.tims=startopt.tims-keepopt.tims;
end

Sprime=(Splus-Smin)/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'rep')
	GTGM=INV.ATA.*SST+opt.mu^2*INV.REGREP;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REGREP*startopt.tims;	
else
	GTGM=INV.ATA.*SST+opt.mu^2*INV.REG;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REG*startopt.tims;
end
testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	deltau=pinv(GTGML)*righths;
	newtime=startopt.tims+deltau; 
	if strcmp(startopt.type,'dep')
		newtime(newtime<0)=0;
	elseif strcmp(startopt.type,'rep')
		newtime(newtime-keepopt.tims < 54)=keepopt.tims(newtime-keepopt.tims < 54) + 54;
	elseif strcmp(startopt.type,'apd')
		newtime(newtime<54)=54;
	end	
	
	testopt.tims=newtime;
	if strcmp(startopt.type,'apd')
		testopt.tims=testopt.tims+keepopt.tims;
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);
	else
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);		
	end
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			opt.tims=testopt.tims; 
			score=1;
		end
		break;
	end
end  % end of inner loop
opt.lambopt=lamb/4;

%% =======================================================================
function TST=gettres_pot(INV,opt,depopt,repopt,otheropt)

global lpass
if strcmp(opt.type,'notch')
	S=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,opt.pot,otheropt.pot,INV.mode);
else %type == ' ampl',
	S=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,otheropt.pot,opt.pot,INV.mode);
end
TST.PHIA=lowpassma(INV.AMA*S(:,opt.usetimes),lpass);
% TST.PHIA=baselinecor(lowpassma(INV.AMA*S(:,opt.usetimes),lpass));
TST.RES=INV.PHIREF(:,opt.usetimes)-TST.PHIA;
TST.rd=norm(TST.RES,'fro')/norm(INV.PHIREF(:,opt.usetimes),'fro'); %NOTE: unfiltered rd
TST.reg=norm(INV.REGOP*opt.pot)/1000;
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;

%% =======================================================================
function TST=gettres_v(INV,opt,keepopt,notchopt,amplopt)

global lpass
if strcmp(opt.type,'dep')
	TST.S=getS(INV.T,opt.tims,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
else %pol=='rep',
	TST.S=getS(INV.T,keepopt.tims,opt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
end
% TST.PHIA=baselinecor(lowpassma(INV.AMA*TST.S,lpass));
TST.PHIA=lowpassma(INV.AMA*TST.S,lpass);
TST.RES=INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));
TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
if strcmp(opt.type,'rep')
	TST.reg=norm(INV.REGOPREP*opt.tims)/1000;
else
	TST.reg=norm(INV.REGOP*opt.tims)/1000;
end
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);

%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM,useSurfLapl)

if useSurfLapl==1
	if max(max(abs(GEOM.VER)))>1
		REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
	else
		REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
	end
	REGOPREP=REGOP;
% 
%     REGOPREP=calcAnisADJ(GEOM,0.25);
% 	REGOPREP(GEOM.ADJsurf>0)=0;
% 	for i=1:length(REGOPREP)
% 		a=[(1:length(REGOPREP))' REGOPREP(:,i)];a=sortrows(a,2);
% 		a(a(:,2)==0,:)=[];a=a(8:end,:);
% 		REGOPREP(a(:,1),i)=0;		REGOPREP(i,a(:,1))=0;
% 	end
% 	REGOPREP(GEOM.ADJsurf>0)=GEOM.ADJsurf(GEOM.ADJsurf>0);
% % 	REGOPREP=graphdist(REGOPREP);
% 	if max(max(abs(GEOM.VER)))>1
% 		REGOPREP=REGOPREP/1000;
% 	end
% 	REGOPREP=REGOPREP.^(-2);
% 	REGOPREP(REGOPREP==Inf)=0;
% 	scal=sum(REGOPREP,2);
% 	scal=1./scal;
% 	REGOPREP=diag(scal)*REGOPREP;
% 	REGOPREP=(REGOPREP-eye(size(REGOPREP)));% scales the result to be of the same order of SL	
% 	fact=[];
% 	for i=1:length(REGOPREP)
% 		fact(i)=REGOP(i,i)/REGOPREP(i,i);
% 	end
% 	REGOPREP=REGOPREP*mean(fact);
%     % 59.7
% 	REGOPREP=597*(REGOPREP-eye(size(REGOPREP)));% scales the result to be of the same order of SL
% REGOP(1,REGOP(1,:)~=0)
% REGOPREP(1,REGOPREP(1,:)~=0)
% REGOPREP=REGOP;
else %if useSurfLapl==2
	disp('inverse 1.distance^2 regularization');
% 	REGOP=GEOM.DIST2W;
%     REGOP=GEOM.ANISDIST;
%     ADJ=GEOM.ADJ; ADJ(ADJ>50)=0;
%     REGOP=calcAnisADJ(GEOM.VER,GEOM.ITRI,ADJ,GEOM.RVER,GEOM.LVER,2.5);
    REGOP=calcAnisADJ(GEOM,2.5);
	if max(max(abs(GEOM.VER)))>1
		REGOP=REGOP/1000;
	end
	REGOP=REGOP.^(-2);
	REGOP(REGOP==Inf)=0;
	scal=sum(REGOP,2);
	scal=1./scal;
	REGOP=diag(scal)*REGOP;
    % 59.7
	REGOP=70*(REGOP-eye(size(REGOP)));% scales the result to be of the same order of SL
end

%% =======================================================================
   
function S=getS(T,dep,rep,p,notchpot,scaleAmpl,mode)
global dodeprep
if dodeprep
    if size(T,2) < max(rep) + 50
        Tt=ones(length(dep),1)*([0:max(rep) + 50]);     
         S=getSmode(Tt,dep,rep,p,scaleAmpl,4);
         S=S(:,1:size(T,2));
    else
        S=getSmode(T,dep,rep,p,scaleAmpl,4);
    end
else
% S=getSvplateau(T,dep,rep,p,notchpot,mode);
S=getSmode(T,dep,rep,p,notchpot,mode);
end