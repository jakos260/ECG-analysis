function PHI=calcEGM(varargin)

%====================================================================
%	P.M. van Dam
%	Viatron Systems Incorporated
%
%====================================================================

if length(varargin)<3
	error('Not enough parameters!!!');
else
	S=varargin{1};
	obsP=varargin{2};
	type=varargin{3};
	doplot=0;
	doCavities=1;
	if length(varargin)>3
		doplot=varargin{4};
	end
	if length(varargin)>4
		doCavities=varargin{5};
	end

end

global ATRIA
global VENTR
global TORSO

sigmaG = 1;		% general
sigmaL = .2;	% lungs
sigmaB = 3;		% blood
if isfield(ATRIA,'sigmaE') % electrode
    sigmaE = ATRIA.sigmaE;
else
    sigmaE =1;
end
sigmaK = sigmaG;

% Determine where the observation point is
pa= find(ATRIA.VER(:,1)  ==obsP(1) & ATRIA.VER(:,2)  ==obsP(2) & ATRIA.VER(:,3)  ==obsP(3));
pv= find(VENTR.VER(:,1)  ==obsP(1) & VENTR.VER(:,2)  ==obsP(2) & VENTR.VER(:,3)  ==obsP(3));
plc=find(ATRIA.LVER(:,1) ==obsP(1) & ATRIA.LVER(:,2) ==obsP(2) & ATRIA.LVER(:,3) ==obsP(3));
prc=find(ATRIA.RVER(:,1) ==obsP(1) & ATRIA.RVER(:,2) ==obsP(2) & ATRIA.RVER(:,3) ==obsP(3));
prl=find(TORSO.RLVER(:,1)==obsP(1) & TORSO.RLVER(:,2)==obsP(2) & TORSO.RLVER(:,3)==obsP(3));
pll=find(TORSO.LLVER(:,1)==obsP(1) & TORSO.LLVER(:,2)==obsP(2) & TORSO.LLVER(:,3)==obsP(3));
pt=find(TORSO.TVER(:,1)  ==obsP(1) & TORSO.TVER(:,2) ==obsP(2) & TORSO.TVER(:,3) ==obsP(3));
pel=[];

if type == 1
	if isfield(ATRIA,'AMA_E')
		pel=find(ATRIA.EVER(:,1)==obsP(1) & ATRIA.EVER(:,2)==obsP(2) & ATRIA.EVER(:,3)==obsP(3));
	end
	if ~isempty(pa) 
		PHI=ATRIA.AMA_A*S;		PHI=PHI(pa,:);	return
	elseif ~isempty(plc)
		PHI=ATRIA.AMA_LC*S;		PHI=PHI(plc,:);	return
	elseif ~isempty(prc)
		PHI=ATRIA.AMA_RC*S;		PHI=PHI(prc,:);	return
	elseif ~isempty(prl)
		PHI=TORSO.AMA_A_RL*S;	PHI=PHI(prl,:); return
	elseif ~isempty(pll)
		PHI=TORSO.AMA_A_LL*S;	PHI=PHI(pll,:);	return
	elseif ~isempty(pel)
		PHI=ATRIA.AMA_E*S;		PHI=PHI(pel,:);	return
	elseif ~isempty(pt)
		PHI=TORSO.AMA_A*S;		PHI=PHI(pt,:);	return
	elseif ~isempty(pv)
		PHI=VENTR.AMA_A*S;		PHI=PHI(pv,:);	return
	end
else
	if isfield(VENTR,'AMA_E')
		pel=find(ATRIA.EVER(:,1)==obsP(1) & ATRIA.EVER(:,2)==obsP(2) & ATRIA.EVER(:,3)==obsP(3));
	end
	if ~isempty(pv)
		PHI=VENTR.AMA_V*S;		PHI=PHI(pv,:);	return
	elseif ~isempty(plc)
		PHI=VENTR.AMA_LC*S;		PHI=PHI(plc,:);	return
	elseif ~isempty(prc)
		PHI=VENTR.AMA_RC*S;		PHI=PHI(prc,:);	return
	elseif ~isempty(prl)
		PHI=TORSO.AMA_V_RL*S;	PHI=PHI(prl,:); return
	elseif ~isempty(pll)		
		PHI=TORSO.AMA_V_LL*S;	PHI=PHI(pll,:);	return
	elseif ~isempty(pel)
		PHI=VENTR.AMA_E*S;		PHI=PHI(pel,:);	return		
	elseif ~isempty(pt)
		PHI=TORSO.AMA_V*S;		PHI=PHI(pt,:);	return
	elseif ~isempty(pa)
		PHI=ATRIA.AMA_V*S;		PHI=PHI(pa,:);	return
	end
end


if sum(solida(ATRIA.VER,ATRIA.ITRI,obsP)) > 0.1 
	% Within the source is undefined
	PHI=100;
	return
elseif sum(solida(VENTR.VER,VENTR.ITRI,obsP)) > 0.1
	% Within the source is undefined
	PHI=200;
	return
elseif sum(solida(TORSO.TVER,TORSO.TITRI,obsP)) < 0.1
	% Outside the body
	PHI=0;
	return;
end

if sum(solida(ATRIA.LVER,ATRIA.LITRI,obsP)) > 0.1 & isfield(ATRIA,'AMA_LC')
	if doCavities==0
		PHI=400;	% Only potentials outside the heart need to be calculated 
		return
	end
	sigmaK=sigmaB;%+sigmaG;
elseif sum(solida(ATRIA.RVER,ATRIA.RITRI,obsP)) > 0.1 & isfield(ATRIA,'AMA_RC')
	if doCavities==0
		PHI=400;	% Only potentials outside the heart need to be calculated 
		return
	end
	sigmaK=sigmaB;%+sigmaG;
elseif sum(solida(TORSO.LLVER,TORSO.LLITRI,obsP)) > 0.1 & isfield(TORSO,'AMA_A_LL')
% 	sigmaK=sigmaL;
	PHI=500;	% Nealy never an electrode is positioned in the lungs
	return
elseif sum(solida(TORSO.RLVER,TORSO.RLITRI,obsP)) > 0.1 & isfield(TORSO,'AMA_A_RL')
	PHI=500;	% Nealy never an electrode is positioned in the lungs
	return
% 	sigmaK=sigmaL;
end

if doplot, figure(1234);clf; hold on; end

if type==1
    PHI=40*sa(ATRIA.VER,ATRIA.ITRI,obsP)*S./sigmaK;		
    if doplot, plot(lowpassma(PHI,3),'b','linewidth',2);lt='atria';end
    if isfield(ATRIA,'AMA_RC')
        phi=((sigmaB-sigmaG)/sigmaK)*sa(ATRIA.RVER,ATRIA.RITRI,obsP)*(ATRIA.AMA_RC*S);
        if doplot, plot(lowpassma(phi,3),'r','linewidth',2);lt=[lt;'Rcave'];end
        PHI=PHI+phi;
    end
    if isfield(ATRIA,'AMA_LC')
        phi=((sigmaB-sigmaG)/sigmaK)*sa(ATRIA.LVER,ATRIA.LITRI,obsP)*(ATRIA.AMA_LC*S);
        if doplot, plot(lowpassma(phi,3),'r:','linewidth',2);lt=[lt;'Lcave'];end
        PHI=PHI+phi;
	 end
    if isfield(TORSO,'AMA_A')
        phi=((sigmaG)/sigmaK)*sa(TORSO.TVER,TORSO.TITRI,obsP)*(TORSO.AMA_A*S);
        if doplot, plot(lowpassma(phi,3),'g','linewidth',2);lt=[lt;'torso'];end
        PHI=PHI+phi;
	end
    if isfield(TORSO,'AMA_A_RL')
        phi=((sigmaL-sigmaG)/(sigmaK))*sa(TORSO.RLVER,TORSO.RLITRI,obsP)*(TORSO.AMA_A_RL*S);
        if doplot, plot(lowpassma(phi,3),'m','linewidth',2);lt=[lt;'Rlung'];end
        PHI=PHI+phi;
    end
    if isfield(TORSO,'AMA_A_LL')
        phi=((sigmaL-sigmaG)/(sigmaK))*sa(TORSO.LLVER,TORSO.LLITRI,obsP)*(TORSO.AMA_A_LL*S);
        if doplot, plot(lowpassma(phi,3),'m:','linewidth',2);lt=[lt;'Llung'];end
        PHI=PHI+phi;
	 end 
    if isfield(ATRIA,'AMA_E')
        phi=((sigmaE-sigmaG)/sigmaK)*sa(ATRIA.EVER,ATRIA.EITRI,obsP)*(ATRIA.AMA_E*S);
        if doplot, plot(lowpassma(phi,3),'y','linewidth',2);lt=[lt;'elect'];end
        PHI=PHI+phi;
	 end	 
    if doplot, plot(lowpassma(PHI,4),'k','linewidth',2); lt=[lt;'final']; legend(lt);drawnow; end
elseif type==2
    PHI=40*sa(VENTR.VER,VENTR.ITRI,obsP)*S/sigmaK;
    if doplot, plot(lowpassma(PHI,3),'b');end
    if isfield(VENTR,'AMA_RC')
        phi=((sigmaB-sigmaG)/sigmaK)*sa(ATRIA.RVER,ATRIA.RITRI,obsP)*(VENTR.AMA_RC*S);
        if doplot, plot(lowpassma(phi,3),'r');end
        PHI=PHI+phi;
    end
    if isfield(VENTR,'AMA_LC')
        phi=((sigmaB-sigmaG)/sigmaK)*sa(ATRIA.LVER,ATRIA.LITRI,obsP)*(VENTR.AMA_LC*S);
        if doplot, plot(lowpassma(phi,3),'g');end
        PHI=PHI+phi;
	 end
    if isfield(TORSO,'AMA_V')
        phi=(sigmaG/sigmaK)*sa(TORSO.TVER,TORSO.TITRI,obsP)*(TORSO.AMA_V*S);
        if doplot, plot(lowpassma(phi,3),'m');end
        PHI=PHI+phi;
	 end
    if isfield(TORSO,'AMA_V_RL')
        phi=((sigmaL-sigmaG)/sigmaK)*sa(TORSO.RLVER,TORSO.RLITRI,obsP)*(TORSO.AMA_V_RL*S);
        if doplot, plot(lowpassma(phi,3),'m');end
        PHI=PHI+phi;
    end
    if isfield(TORSO,'AMA_V_LL')
        phi=((sigmaL-sigmaG)/sigmaK)*sa(TORSO.LLVER,TORSO.LLITRI,obsP)*(TORSO.AMA_V_LL*S);
        if doplot, plot(lowpassma(phi,3),'c');end
        PHI=PHI+phi;
	end 
    if isfield(VENTR,'AMA_E')
        phi=((sigmaE-sigmaG)/sigmaK)*sa(ATRIA.EVER,ATRIA.EITRI,obsP)*(VENTR.AMA_E*S);
        if doplot, plot(lowpassma(phi,3),'y','linewidth',2);end
        PHI=PHI+phi;
	 end	 
    if doplot, plot(lowpassma(PHI,4),'k','linewidth',2); drawnow; end
end

%************************************************************

function rowb=sa(VER,ITRI,obsP)

[OMEGA,index]=dsa(VER,ITRI,obsP,.1);
rowb=zeros(1,length(VER));
for j=1:length(ITRI),
    rowb(ITRI(j,:))=rowb(ITRI(j,:)) + OMEGA(:,j)';
end
rowb=rowb/(4*pi);

%************************************************************
