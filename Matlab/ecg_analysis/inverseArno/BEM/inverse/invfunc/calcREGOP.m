%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM,useSurfLapl)

% if isfield(GEOM,'LAPL')
%     REGOP2 = 2.0 * GEOM.LAPL;
%     REGOPREP2 = REGOP;
% end
if useSurfLapl==1
	if max(max(abs(GEOM.VER)))>1
		REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
	else
		REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
    end
%     if isfield(GEOM,'LAPL')
%         REGOP = REGOP - GEOM.LAPL;
%     end
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
	% scal=sum(REGOPREP,2);
	% scal=1./scal;
	% REGOPREP=diag(scal)*REGOPREP;
	% REGOPREP=(REGOPREP-eye(size(REGOPREP)));% scales the result to be of the same order of SL	
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
    [a,ORDER] = graphdist(GEOM.ITRI);
 	REGOP=GEOM.DIST2W;
%     REGOP=GEOM.ANISDIST;
%     ADJ=GEOM.ADJ; ADJ(ADJ>50)=0;
%     REGOP=calcAnisADJ(GEOM.VER,GEOM.ITRI,ADJ,GEOM.RVER,GEOM.LVER,2.5);
%     REGOP=calcAnisADJ(GEOM,2.5);
	if max(max(abs(GEOM.VER)))>1
		REGOP=3*REGOP/1000;
	end
	REGOP=REGOP.^(-1);
    REGOPREP = REGOP;
    if max(max(abs(GEOM.VER)))>1
		B=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
	else
		B=surflapl(GEOM.VER,GEOM.ITRI,1);
    end

    REGOP(GEOM.DIST>10)=0;
    REGOP(ORDER<=4)=0;
    B2=diag(-0.3*sum(REGOP));
    REGOP(B~=0)=B(B~=0);
    REGOP =  REGOP + B2;
%     REGOPREP(GEOM.DIST>20)=0;
%     REGOPREP(ORDER<=1)=0;
%     REGOPREP(B~=0)=B(B~=0);
    REGOPREP=B;

% 	REGOP(REGOP==Inf)=0;
% 	scal=sum(REGOP,2);
% 	scal=1./scal;
% 	REGOP=diag(scal)*REGOP;
%     % 59.7
% 	REGOP=70*(REGOP-eye(size(REGOP)));% scales the result to be of the same order of SL
end
% --- normalizacja macierzy Laplasjanu ---
% 1. Pobieramy wartości z przekątnej
diag_vals = diag(REGOPREP);
% 2. Zabezpieczenie przed ewentualnym zerem (izolowane węzły na siatce)
diag_vals(diag_vals == 0) = 1; 
% 3. Skalujemy wiersze przez odwrotność wartości bezwzględnej z przekątnej
scal_rep = 1 ./ abs(diag_vals);
REGOPREP = diag(scal_rep) * REGOPREP;

% To samo dla macierzy depolaryzacji
diag_vals_dep = diag(REGOP);
diag_vals_dep(diag_vals_dep == 0) = 1;
scal_dep = 1 ./ abs(diag_vals_dep);
REGOP = diag(scal_dep) * REGOP;
end