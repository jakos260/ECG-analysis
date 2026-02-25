function varargout = computeEGM(varargin)


if length(varargin)<2
    error('Not enough parameters!!!');
else
    GEOM = varargin{1};
    obsP=varargin{2};
    doplot=1;
    S=[];
    if length(varargin)>2
        S=varargin{3};
    end
    if length(varargin)>3
        doplot=varargin{4};
    end
end

sigmaG = 1;		% general
sigmaL = .2;	% lungs
sigmaB = 3;		% blood
sigmaE = 100;
sigmaK = sigmaG;

% Determine where the observation point is
pa = find(GEOM.VER(:,1)  ==obsP(1) & GEOM.VER(:,2)   ==obsP(2) & GEOM.VER(:,3)  ==obsP(3));
plc= find(GEOM.LVER(:,1) ==obsP(1) & GEOM.LVER(:,2)  ==obsP(2) & GEOM.LVER(:,3) ==obsP(3));
prc= find(GEOM.RVER(:,1) ==obsP(1) & GEOM.RVER(:,2)  ==obsP(2) & GEOM.RVER(:,3) ==obsP(3));
prl= find(GEOM.RLVER(:,1)==obsP(1) & GEOM.RLVER(:,2) ==obsP(2) & GEOM.RLVER(:,3)==obsP(3));
pll= find(GEOM.LLVER(:,1)==obsP(1) & GEOM.LLVER(:,2) ==obsP(2) & GEOM.LLVER(:,3)==obsP(3));
pt = find(GEOM.TVER(:,1)  ==obsP(1) & GEOM.TVER(:,2) ==obsP(2) & GEOM.TVER(:,3) ==obsP(3));
pel=[];

AMA =[];
% 	if isfield(ATRIA,'AMA_E')
% 		pel=find(ATRIA.EVER(:,1)==obsP(1) & ATRIA.EVER(:,2)==obsP(2) & ATRIA.EVER(:,3)==obsP(3));
% 	end
if ~isempty(pa)
    AMA = GEOM.AMAH;        if ~isempty(S), PHI = AMA(pa,:) * S;  end;
elseif ~isempty(plc)
    AMA = GEOM.AMALCAV;     if ~isempty(S), PHI = AMA(plc,:) * S;  end;
elseif ~isempty(prc)
    AMA=GEOM.AMARCAV;		if ~isempty(S), PHI = AMA(prc,:) * S;  end;
elseif ~isempty(prl)
    AMA=GEOM.AMARL*S;       if ~isempty(S), PHI = AMA(prl,:) * S;  end;
elseif ~isempty(pll)
    AMA=GEOM.AMALL;         if ~isempty(S), PHI = AMA(pll,:) * S;  end;
    % 	elseif ~isempty(pel)
    % 		AMA=GEOM.AMA_E*S;		PHI = AMA(pel,:) *S;	return
elseif ~isempty(pt)
    AMA=GEOM.AMA*S;         if ~isempty(S), PHI = AMA(pt,:) * S;  end;
end
if exist('PHI')
    varargout{1}  = AMA;
    if nargout ==2
        varargout{2} = PHI;
    end
end

if sum(solida(GEOM.VER,GEOM.ITRI,obsP)) > 0.1
    % Within the source is undefined
    AMA = [];
    PHI = [];
    return
elseif sum(solida(GEOM.TVER,GEOM.TITRI,obsP)) < 0.1
    % Outside the body
    AMA = [];
    PHI=[];
    return;
end


if doplot, figure(1234);clf; hold on; end

AMA = 40*sa(GEOM.VER,GEOM.ITRI,obsP)./sigmaK;
if ~isempty(S)
    PHI = AMA * S;
    if doplot, plot(lowpassma(PHI,3),'b','linewidth',2);lt='heart';end
end
A = AMA;

if isfield(GEOM,'AMARCAV')
    AMARCAV=((sigmaB-sigmaG)/sigmaK)*sa(GEOM.RVER,GEOM.RITRI,obsP)*(GEOM.AMARCAV);
    A = A + AMARCAV;
    if ~isempty(S)
        
        PHIRCAV = AMARCAV * S;
        if doplot, plot(lowpassma(PHIRCAV,3),'r','linewidth',2);lt=[lt;'Rcave'];end
        PHI=PHI + PHIRCAV;
    end
end
if isfield(GEOM,'AMALCAV')
    AMALCAV=((sigmaB-sigmaG)/sigmaK)*sa(GEOM.LVER,GEOM.LITRI,obsP)*(GEOM.AMALCAV);
    A = A + AMALCAV;
    if ~isempty(S)
        
        PHILCAV = AMALCAV * S;
        if doplot, plot(lowpassma(PHIRCAV,3),'r:','linewidth',2);lt=[lt;'Lcave'];end
        PHI=PHI + PHILCAV;
    end
end
if isfield(GEOM,'AMA')
    AMAT=((sigmaG)/sigmaK) * sa(GEOM.TVER,GEOM.TITRI,obsP)*GEOM.AMAORG;
    A = A + AMAT;
    if ~isempty(S)
        
        PHIT = AMAT * S;
        if doplot, plot(lowpassma(PHIT,3),'g','linewidth',2);lt=[lt;'torso'];end
        PHI = PHI + PHIT;
    end
end
if isfield(GEOM,'AMARL')
    AMARL=((sigmaL-sigmaG)/(sigmaK))*sa(GEOM.RLVER,GEOM.RLITRI,obsP)*(GEOM.AMARL);
    A = A + AMARL;
    if ~isempty(S)
        
        PHIRL = AMARL * S;
        if doplot, plot(lowpassma(PHIRL,3),'m','linewidth',2);lt=[lt;'Rlung'];end
        PHI = PHI + PHIRL;
    end
end
if isfield(GEOM,'AMALL')
    AMALL=((sigmaL-sigmaG)/(sigmaK))*sa(GEOM.LLVER,GEOM.LLITRI,obsP)*(GEOM.AMALL);
    A = A + AMALL;
    if ~isempty(S)        
        PHILL = AMALL * S;
        if doplot, plot(lowpassma(PHIRL,3),'m:','linewidth',2);lt=[lt;'Rlung'];end
        PHI = PHI + PHILL;
    end
end
%     if isfield(ATRIA,'AMA_E')
%         phi=((sigmaE-sigmaG)/sigmaK)*sa(ATRIA.EVER,ATRIA.EITRI,obsP)*(ATRIA.AMA_E*S);
%         if doplot, plot(lowpassma(phi,3),'y','linewidth',2);lt=[lt;'elect'];end
%         PHI = PHI + PHILL;
%     end
if ~isempty(S) && doplot, plot(lowpassma(PHI,4),'k','linewidth',2); lt=[lt;'final']; legend(lt);drawnow; end

varargout{1}  = A;
if nargout ==2
    varargout{2} = PHI;
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
