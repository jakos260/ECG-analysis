% function S = getSmode(T,dep,rep,SPECS,mode,varargin)
%
% This function gives the simulated Transmembrane potentials [TMP's] for the
% given depolarization & repolarization times, with or without the
% Anti-Hakkel --> if varargin = GEOM, Anti-Hakkel is used.
%
% INPUT:
% T     = matrix of # vertices * timescale QRT
% dep   = depolarization times
% rep   = repolarization times
% SPECS = Contains information about slopes in RMS of signal. Needed to construct TMP's.
% mode  = use only depolarization or also repolarization
%
% OUTPUT:
% S     = Simulated TMP's.
%
% variant of getS dedicated to ventricles function:
% [S = getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)]; PM van Dam; 2008-10-11
%
% Version 1: 01-apr-2015

% New version Jeanne van der Waal; 2018-16-5
% new mode: don't use GEOM.SPECS parameters for repolarization


function S = getSmodeJ(T,dep,rep,SPECS,mode,varargin)
new_mode=0;     % default the 'old' mode
if ~isfield(SPECS,'useCumsum')
    SPECS.useCumsum=1;
end
% check the dimensions of dep and rep:
nt = size(T,2);
if size(dep,1) < size(dep,2), dep = dep'; end
if size(rep,1) < size(rep,2), rep = rep'; end

%% Anti-Hakkel: if 6 inputs are given instead of 5 [GEOM]
if length(varargin) == 1
    if length(varargin{1})>1
        neigh   = varargin{1};
        S       = getSpc(T,dep,rep,SPECS,mode,neigh);               % The adjacency matrix with graphdist(ITRI) [GEOM.neigh] is used.
        return;
    else
        new_mode=varargin{1};
    end
    
end

%% NO Anti-Hakkel:

% Depolarization:
TDEP = T - dep*ones(1,nt);
if size(SPECS.depSlope,1) == 1
    S = 1./(1+exp(-SPECS.depSlope.*TDEP));                      % describe TMP's as sigmoid with a specified depolarization slope.
    Sdep=S;
else
    S = 1./(1+exp(bsxfun(@times,TDEP,-SPECS.depSlope(:))));     % @times --> Array multiply: multiple depolarization slopes for TMP's.
end

% Repolarization:
if mode >= 4
    warning('off','MATLAB:interp1:UsePCHIP')
    if new_mode==0
        TREP = T-rep*ones(1,nt)-SPECS.repCorrection;
        
        if SPECS.useCumsum
            Y = (SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP));
            Y = 1 - bsxfun(@times,cumsum(Y,2),1./sum(Y,2));
            
        else
            Y = -(SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP));
            Y = bsxfun(@times,Y,1./Y(:,1));
        end
        
        % apply the up-slope: logistic shape;
        S = S .* Y;
        if ~isempty(SPECS.scaleAmpl)
            S = bsxfun(@rdivide,S,max(S,[],2)./SPECS.scaleAmpl);
        else
            S = bsxfun(@rdivide,S,max(S,[],2));
        end
    elseif new_mode==1                          % added mode - JvdW for TMP in repolarization
        % mode 1 = platau until a third of amp (so adapt steepness of plateauslope with time)
        TREP = T-rep*ones(1,nt);     % RepCorrection?? For now - take 10 ms
        
        P3=(1./(1+exp(-0.02.*TREP)))./(1+exp(0.09.*TREP));          % Phase 3 of action potential
        P3=1-(cumsum(P3,2).*(1./sum(P3,2)));
        P3(P3>0.669)=1;
        P23=zeros(size(P3));%P2_2=P2;
        for x=1:size(TREP,1)
            deptime=round(dep(x));
            endplateau=find(P3(x,:)<0.67);
            endplateau=endplateau(1);
            linslope=resample([1 0.67],1,2,endplateau-deptime-1);
            %             step=0.33/(endplateau-deptime);
            %             endtime=size(T,2)-endplateau;
            %             P2(x,:)=[ones(1,deptime) linslope' ones(1,endtime+1)];    %
            P23(x,:)=[ones(1,deptime) linslope' P3(x,endplateau:end)];% phase 2 of action potential +phase 3
            % For P2: try quadratic function or exponential function
            % instead of linear decrease:
            %             plat=1/1.5^(deptime/(deptime-endplateau))*((1.5^(1/(deptime-endplateau))).^(deptime:endplateau-1));
            %             P2_2(x,:)=[ones(1,deptime) plat ones(1,endtime)];    % phase 2 of action potential with exp func
        end
        %         P23=P2.*P3;
        S=Sdep.*P23;
        %         S = bsxfun(@rdivide,S,max(S,[],2));
        %         S_2=Sdep.*P23_2;
    elseif new_mode==2                          % mode 2 = 'rescaling' of entire curve
        
        %         TREP = T-rep*ones(1,nt)-15;     % RepCorrection?? For now - take 10 ms
        
        for x=1:size(dep,1)
            deptime=round(dep(x));
            reptime=round(rep(x));
            AP=[0 0.5 1 0.998:-0.002:0.994 0.152*(-0.001*((-9.5:0.04:-1.5).^(3))+1)+0.71 0.01*(-1*(0.5:0.005:1).^3+2)+0.8437 0.8536:-0.1/420:0.75 0.085*(-1*(0.5:0.004:2.08).^3+2)+0.59];
            AP(AP<0)=0;
            rt=1; ad=0;
            while rt                    % Zoek juiste lengte om te resamplen --> hoogte van AP moet 0.2 zijn op reptime
                lengthAP=reptime-deptime+ad;
                ad=ad+1;
                APtemp=resample(AP,1,length(AP),lengthAP);
                if APtemp(reptime-deptime)>0.19
                    rt=0;
                    APnew=APtemp;
                end
            end
            if length(APnew)+deptime>nt
                stemp=[zeros(1,deptime) APnew'];
                stemp(nt+1:end)=[];
                S(x,:)=stemp;
            else
                S(x,:)=[zeros(1,deptime) APnew' zeros(1,nt-length(APnew)-deptime)];
            end
        end
    elseif new_mode==3          % mode 3 = 'rescaling' of entire curve but WITHOUT initial (phase 1) dip
        % modeled with parabola
        for x=1:size(dep,1)
            deptime=round(dep(x));
            reptime=round(rep(x));
            AP=[0 0.5 ones(1,10) -((1:800)/800).^2+1];
            rt=1; ad=0;
            while rt                    % Zoek juiste lengte om te resamplen --> hoogte van AP moet 0.2 zijn op reptime
                lengthAP=reptime-deptime+ad;
                ad=ad+1;
                APtemp=resample(AP,1,length(AP),lengthAP);
                if APtemp(reptime-deptime)>0.19
                    rt=0;
                    APnew=APtemp;
                end
            end
            if length(APnew)+deptime>nt
                stemp=[zeros(1,deptime) APnew'];
                stemp(nt+1:end)=[];
                S(x,:)=stemp;
            else
                S(x,:)=[zeros(1,deptime) APnew' zeros(1,nt-length(APnew)-deptime)];
                S(S<0)=0;
            end
            
        end
    elseif new_mode==4          % mode 4 = lengthening/shortening of only plateauphase, but
        % keeping slope of rep-phase relatively
        % equally long
        for x=1:size(dep,1)
            deptime=round(dep(x));
            reptime=round(rep(x));
            apd=reptime-deptime;
            x1=2.1:15/(apd*1.05)*3:15;
            y=(-x1.^3+15^3)/(15^3+20);
            step=(y(1)-1)/(2*apd/3);
            x2=1:step:y(1)-step;
            AP=[zeros(1,deptime) x2 y zeros(1,nt-length([y x2])-deptime)];
            S(x,:)=AP;
        end
    elseif new_mode==5          % mode 5 = reptime is always around 14% of AP-amplitude
                                % modeled with x^3 formula
        for x=1:size(dep,1)
            deptime=round(dep(x));
            reptime=round(rep(x));
            apd=reptime-deptime;
            x1=0.1:15/(apd*1.05):15;
            y=(-x1.^3+15^3)/(15^3+20);
            AP=[zeros(1,deptime) y zeros(1,nt-length(y)-deptime)];
            AP=AP(1:nt);
            S(x,:)=AP;
            S = bsxfun(@rdivide,S,max(S,[],2));
        end
    elseif new_mode==6          % mode = 6
        AP=[0 0.5 1 1/0.8625*([0.01*(-1*(0.5:0.005:1).^3+2)+0.8437 0.8536:-0.1/420:0.75 0.085*(-1*(0.5:0.004:2.08).^3+2)+0.59])];
        for x=1:size(dep,1);
            deptime=round(dep(x));
            reptime=round(rep(x));
            rt=1; ad=0;
            while rt                    % Zoek juiste lengte om te resamplen --> hoogte van AP moet 0.2 zijn op reptime
                lengthAP=reptime-deptime+ad;
                ad=ad+1;
                APtemp=resample(AP,1,length(AP),lengthAP);
                if APtemp(reptime-deptime)>0.19
                    rt=0;
                    APnew=APtemp;
                end
            end
            if length(APnew)+deptime>nt
                stemp=[zeros(1,deptime) APnew'];
                stemp(nt+1:end)=[];
                S(x,:)=stemp;
                S(S<0)=0;
            else
                S(x,:)=[zeros(1,deptime) APnew' zeros(1,nt-length(APnew)-deptime)];
                S(S<0)=0;
            end
        end
    end
end