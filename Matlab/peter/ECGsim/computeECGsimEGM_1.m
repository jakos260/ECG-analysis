function varargout = computeECGsimEGM(varargin)

sigmaG = 1;		% general
sigmaL = .2;	% lungs

sigmaB = 3;		% blood
sigmaF = 0.16;  % fat
sigmaR = 0.06;  % ribcage
sigmaH = 1;     % heart= source


if nargin<4
    error('Not enough parameters!!!');
else
    GEOM = varargin{1};
    VOL = varargin{2};
    obsP=varargin{3};
    doplot=0;
    S=varargin{4};
    SAH = S;
    if length(varargin)>4
        doplot=varargin{5};
    end
    pp=5;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'fat'
                    sigmaF=varargin{pp+1};pp=pp+2;
                case 'ribcage'
                    sigmaR=varargin{pp+1};pp=pp+2;
                case 'blood'
                    sigmaB=varargin{pp+1};pp=pp+2;
                case 'lungs'
                    sigmaL=varargin{pp+1};pp=pp+2;
                case 'plot'
                    doplot=varargin{pp+1};pp=pp+2;
                case 'antihak'
                    SAH = varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        else
            error('unknown parameter');
        end
    end
end




small = 1e-5;
[aProjPnt,aTri,aLambda,aMu]=pointOnGeometry(    GEOM.atria.VER, GEOM.atria.ITRI,obsP);
if norm(aProjPnt-obsP) < small && ~isempty(VOL.ATRIA)
    itri= GEOM.atria.ITRI(aTri,:);
    cavityStr='on atria';
    AMA = VOL.ATRIA(itri(1),:) + (VOL.ATRIA(itri(2),:) -VOL.ATRIA(itri(1),:)) * aLambda +...
        (VOL.ATRIA(itri(3),:) -VOL.ATRIA(itri(1),:)) * aMu;
else
    [vProjPnt,vTri,vLambda,vMu]=pointOnGeometry(    GEOM.ventr.VER, GEOM.ventr.ITRI,obsP);
    if norm(vProjPnt-obsP) < small
        itri= GEOM.ventr.ITRI(vTri,:);
        cavityStr='on ventricle';
        AMA = VOL.VENTRICLES(itri(1),:) + (VOL.VENTRICLES(itri(2),:) -VOL.VENTRICLES(itri(1),:)) * vLambda +...
            (VOL.VENTRICLES(itri(3),:) -VOL.VENTRICLES(itri(1),:)) * vMu;
    else
        [lcProjPnt,lcTri,lcLambda,lcMu]=pointOnGeometry(GEOM.lcav.VER,  GEOM.lcav.ITRI,obsP);
        if norm(lcProjPnt-obsP) < small
            itri= GEOM.lcav.ITRI(lcTri,:);
            cavityStr='on left cavity';
            AMA = VOL.LCAV(itri(1),:) + (VOL.LCAV(itri(2),:) -VOL.LCAV(itri(1),:)) * lcLambda +...
                (VOL.LCAV(itri(3),:) -VOL.LCAV(itri(1),:)) * lcMu;
        else
            [rcProjPnt,rcTri,rcLambda,rcMu]=pointOnGeometry(GEOM.rcav.VER,  GEOM.rcav.ITRI,obsP);
            if norm(rcProjPnt-obsP) < small
                cavityStr='on right cavity';
                itri= GEOM.rcav.ITRI(rcTri,:);
                AMA = VOL.RCAV(itri(1),:) + (VOL.RCAV(itri(2),:) -VOL.RCAV(itri(1),:)) * rcLambda +...
                    (VOL.RCAV(itri(3),:) -VOL.RCAV(itri(1),:)) * rcMu;
            else
                [rlProjPnt,rlTri,rlLambda,rlMu]=pointOnGeometry(GEOM.rlung.VER, GEOM.rlung.ITRI,obsP);
                if norm(rlProjPnt-obsP) < small
                    cavityStr='on right lung';
                    itri= GEOM.rlung.ITRI(rlTri,:);
                    AMA = VOL.RLUNG(itri(1),:) + (VOL.RLUNG(itri(2),:) -VOL.RLUNG(itri(1),:)) * rlLambda +...
                        (VOL.RLUNG(itri(3),:) -VOL.RLUNG(itri(1),:)) * rlMu;
                else
                    [llProjPnt,llTri,llLambda,llMu]=pointOnGeometry(GEOM.llung.VER, GEOM.llung.ITRI,obsP);
                    if norm(llProjPnt-obsP) < small
                        cavityStr='on left lung';
                        itri= GEOM.llung.ITRI(llTri,:);
                        AMA = VOL.LLUNG(itri(1),:) + (VOL.LLUNG(itri(2),:) -VOL.LLUNG(itri(1),:)) * llLambda +...
                            (VOL.LLUNG(itri(3),:) -VOL.LLUNG(itri(1),:)) * llMu;
                    else
                        [tProjPnt,tTri,tLambda,tMu]=pointOnGeometry(    GEOM.thorax.VER,GEOM.thorax.ITRI,obsP);
                        if norm(tProjPnt-obsP) < small
                            cavityStr='on torso';
                            itri= GEOM.thorax.ITRI(tTri,:);
                            AMA = VOL.THORAX(itri(1),:) + (VOL.THORAX(itri(2),:) -VOL.THORAX(itri(1),:)) * tLambda +...
                                (VOL.THORAX(itri(3),:) -VOL.THORAX(itri(1),:)) * tMu;
                        else
                            if isfield(GEOM,'ribcage')
                                [ribProjPnt,ribTri,ribLambda,ribMu]=pointOnGeometry(GEOM.ribcage.VER,GEOM.ribcage.ITRI,obsP);
                                if  norm(ribProjPnt-obsP) < small
                                    cavityStr='on ribcage';
                                    itri= GEOM.ribcage.ITRI(ribTri,:);
                                    AMA = VOL.RIBCAGE(itri(1),:) + (VOL.RIBCAGE(itri(2),:) -VOL.RIBCAGE(itri(1),:)) * ribLambda +...
                                        (VOL.RIBCAGE(itri(3),:) -VOL.RIBCAGE(itri(1),:)) * ribMu;
                                end
                            end
                            if isfield(GEOM,'fatpad_1')
                                [fat1ProjPnt,fat1Tri,fat1Lambda,fat1Mu]=pointOnGeometry(GEOM.fatpad_1.VER,GEOM.fatpad_1.ITRI,obsP);
                                if norm(fat1ProjPnt-obsP) < small
                                    cavityStr='on fatpad 1';

                                    itri= GEOM.fatpad_1.ITRI(fat1Tri,:);
                                    AMA = VOL.FATPAD_1(itri(1),:) + (VOL.FATPAD_1(itri(2),:) -VOL.FATPAD_1(itri(1),:)) * fat1Lambda +...
                                        (VOL.FATPAD_1(itri(3),:) -VOL.FATPAD_1(itri(1),:)) * fat1Mu;
                                end
                            end
                            if isfield(GEOM,'fatpad_2')
                                [fat2ProjPnt,fat2Tri,fat2Lambda,fat2Mu]=pointOnGeometry(GEOM.fatpad_2.VER,GEOM.fatpad_2.ITRI,obsP);
                                if norm(fat2ProjPnt-obsP) < small
                                    cavityStr='on fatpad 2';
                                    itri= GEOM.fatpad_2.ITRI(fat2Tri,:);
                                    AMA = VOL.FATPAD_2(itri(1),:) + (VOL.FATPAD_2(itri(2),:) -VOL.FATPAD_2(itri(1),:)) * fat2Lambda +...
                                        (VOL.FATPAD_2(itri(3),:) -VOL.FATPAD_2(itri(1),:)) * fat2Mu;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


if exist('AMA','var')
    disp(['location: ' cavityStr]);
    PHI = lowpassma(AMA * S,10);
    varargout{1} = PHI;
    if nargout > 1
        varargout{2} = cavityStr;
    end
    return;
end

sigmaObs  = sigmaG;


% error('other points are still not implemented correctly')
cavityStr='torso';
insideLL = 0;
insideRL = 0;
insideLC = 0;
insideRC = 0;
insideRB = 0;
insideF1 = 0;
insideF2 = 0;
if sum(solida(VOL.geom.VER,VOL.geom.ITRI,obsP)) > 0.1
    % Within the source is undefined
    PHI = zeros(1,size(S,2));
    varargout{1} = PHI;
    if nargout > 1
        varargout{2} = 'in myocardium';
    end
    return;
elseif sum(solida(GEOM.thorax.VER,GEOM.thorax.ITRI,obsP)) < 0.1
    % Outside the body
    PHI = zeros(1,size(S,2));
    varargout{1} = PHI;
    if nargout > 1
        varargout{2} = 'on thorax';
    end
    return;
elseif sum(solida(GEOM.rlung.VER,GEOM.rlung.ITRI,obsP)) > 1
    sigmaObs  = sigmaL;
    cavityStr='right lung';
    insideRL = 1;
elseif sum(solida(GEOM.llung.VER,GEOM.llung.ITRI,obsP)) > 1
    sigmaObs  = sigmaL;
    cavityStr='left lung';
    insideLL = 1;
elseif sum(solida(GEOM.rcav.VER,GEOM.rcav.ITRI,obsP)) > 1
    sigmaObs  = sigmaB;
    cavityStr='right cavity';
    insideRC = 1;
elseif sum(solida(GEOM.lcav.VER,GEOM.lcav.ITRI,obsP)) > 1
    sigmaObs  = sigmaB;
    cavityStr='left cavity';
    insideLC = 1;
elseif isfield(GEOM,'ribcage') && ...
        sum(solida(GEOM.ribcage.VER,GEOM.ribcage.ITRI,obsP)) > 1
    sigmaObs  = sigmaR;
    cavityStr='ribcage';
    insideRB = 1;
elseif isfield(GEOM,'fatpad_1') && ...
        sum(solida(GEOM.fatpad_1.VER,GEOM.fatpad_1.ITRI,obsP)) > 1
    sigmaObs  = sigmaF;
    cavityStr='fatpad 1';
    insideF1 = 1;
elseif isfield(GEOM,'fatpad_2') && ...
        sum(solida(GEOM.fatpad_2.VER,GEOM.fatpad_2.ITRI,obsP)) > 1
    sigmaObs  = sigmaF;
    cavityStr='fatpad 2';
    insideF2 = 1;
end
disp(['location in : ' cavityStr]);
if doplot, figure(12345);clf; hold on; ncomp=1;end

PHI = 40 * sigmaH * sa(VOL.geom.VER,VOL.geom.ITRI,obsP) * S;
if doplot, plot(lowpassma(PHI / sigmaObs,3),'b','linewidth',1);lt{ncomp}='inf heart';ncomp = ncomp+1;end

%torso
PHIT = (sigmaG-0) * sa(GEOM.thorax.VER,GEOM.thorax.ITRI,obsP) * (VOL.THORAX * SAH);
if doplot, plot(lowpassma(PHIT / sigmaObs,3),'g','linewidth',1);lt{ncomp}= 'torso';ncomp = ncomp+1;end
PHI = PHI + PHIT; % by definition + because always inside

if insideRB
    PHIRIB = (sigmaR-sigmaG) * sa(GEOM.ribcage.VER, GEOM.ribcage.ITRI, obsP ) * ( VOL.RIBCAGE * SAH );
    if doplot, plot(lowpassma(PHIRIB / sigmaObs,3),'c','linewidth',1);lt{ncomp}= 'ribcage';ncomp = ncomp+1;end
    PHI = PHI + PHIRIB;
elseif insideRC
    PHIRCAV = (sigmaB-sigmaG) * (sa(GEOM.rcav.VER,GEOM.rcav.ITRI,obsP) * (VOL.RCAV * SAH));
    if doplot, plot(lowpassma(PHIRCAV / sigmaObs,3),'r','linewidth',1);lt{ncomp}= 'Rcav ';ncomp = ncomp+1;end
    PHI=PHI - PHIRCAV;
elseif insideLC
    PHILCAV = (sigmaB-sigmaG) * sa(GEOM.lcav.VER,GEOM.lcav.ITRI,obsP) * (VOL.LCAV*SAH);
    if doplot, plot(lowpassma(PHILCAV / sigmaObs,3),'r:','linewidth',1);lt{ncomp}= 'Lcav ';ncomp = ncomp+1;end
    PHI=PHI - PHILCAV;
elseif insideRL
    PHIRL= (sigmaL-sigmaG)*sa(GEOM.rlung.VER,GEOM.rlung.ITRI,obsP) * (VOL.RLUNG *SAH);
    if doplot, plot(lowpassma(PHIRL / sigmaObs,3),'m','linewidth',1);lt{ncomp}= 'Rlung';ncomp = ncomp+1;end
    PHI = PHI - PHIRL;
elseif insideLL
    PHILL = (sigmaL-sigmaG)*sa(GEOM.llung.VER,GEOM.llung.ITRI,obsP)*(VOL.LLUNG*SAH); 
    if doplot, plot(lowpassma(PHILL / sigmaObs,3),'m:','linewidth',1);lt{ncomp}= 'Llung';ncomp = ncomp+1;end
    PHI = PHI - PHILL;
elseif insideF1
    PHIFAT1 = (sigmaF-sigmaG) * sa(GEOM.fatpad_1.VER, GEOM.fatpad_1.ITRI, obsP )*(VOL.FATPAD_1*SAH);
    if doplot, plot(lowpassma(PHIFAT1 / sigmaObs,3),'y-','linewidth',1);lt{ncomp}= 'fat_1';ncomp = ncomp+1;end
    PHI = PHI - PHIFAT1;
elseif insideF2
    PHIFAT2 = (sigmaF-sigmaG) * sa(GEOM.fatpad_2.VER, GEOM.fatpad_2.ITRI, obsP )*(VOL.FATPAD_2*SAH);   
    if doplot, plot(lowpassma(PHIFAT2 / sigmaObs,3),'y:','linewidth',1);lt{ncomp}= 'fat_2';ncomp = ncomp+1;end
    PHI = PHI - PHIFAT2;
end

if ~isempty(S) && doplot, plot(lowpassma(PHI / sigmaObs,5),'k','linewidth',1); lt{ncomp}= ['final in ' cavityStr] ; legend(lt);drawnow; end

varargout{1} = lowpassma(PHI / sigmaObs,1);
if nargout > 1
    varargout{2} = cavityStr;
end


%************************************************************

function rowb=sa(VER,ITRI,obsP)

      
OMEGA = distributedSa(VER,ITRI,obsP,0.01);
rowb = zeros(1,length(VER));
for j=1:length(ITRI),
    rowb(ITRI(j,:)) = rowb(ITRI(j,:)) + OMEGA(:,j)';
end
rowb=rowb/(4*pi);
% rowb = rowforw(VER,ITRI,obsP)/2;



%************************************************************

function  SA=distributedSa(VER,ITRI,obs,small)
[nver,~]=size(VER);
VER=VER-ones(nver,1)*obs;

R1=VER(ITRI(:,1),:);
R2=VER(ITRI(:,2),:);
R3=VER(ITRI(:,3),:);
LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];
R2xR3=cross(R2,R3);
ntri=size(ITRI,1);
% blockproducts
block=dots(R1,R2xR3);
DOTS=[dots(R1,R2) dots(R2,R3) dots(R3,R1)];
denom=LR(:,1).*LR(:,2).*LR(:,3) + ...
    LR(:,1).*DOTS(:,2)+LR(:,2).*DOTS(:,3) + LR(:,3).*DOTS(:,1);
sa=-2*atan2(block,denom)';

sa(abs(block)<=eps)=0; % crude approximation
SA  =ones(3,1)*sa/3;
TEST=[1:ntri;sa];
k   =TEST(1,abs(TEST(2,:))>small);

if ~isempty(k)
    % dsa computation,  required for triangles: k
    S1=R2(k,:)-R1(k,:);
    S2=R3(k,:)-R2(k,:);
    S3=R1(k,:)-R3(k,:);
    LS=[norm3d(S1) norm3d(S2) norm3d(S3)];
    DOTS=DOTS(k,:);
    
    R2xR3=R2xR3(k,:);
    R3xR1=cross(R3(k,:),R1(k,:));
    R1xR2=cross(R1(k,:),R2(k,:));
    
    NT    =R1xR2+R2xR3+R3xR1;  % dimension m^2
    asq   =dots(NT,NT);  % dimension m^4
    NT    =(sa(k)'*ones(1,3)).*NT;
    
    % compute the gamma vector
    m=[2 3 1];
    
    NUM = LR(k,:).*LS+DOTS-LR(k,:).*LR(k,:);
    DNOM= LR(k,m).*LS-DOTS+LR(k,m).*LR(k,m);
    
    GAM = ((block(k)*ones(1,3)).*real(log(NUM./DNOM)))./LS;
    
    GV(:,1)=S1(:,1).*GAM(:,1) + S2(:,1).*GAM(:,2) + S3(:,1).*GAM(:,3);
    GV(:,2)=S1(:,2).*GAM(:,1) + S2(:,2).*GAM(:,2) + S3(:,2).*GAM(:,3);
    GV(:,3)=S1(:,3).*GAM(:,1) + S2(:,3).*GAM(:,2) + S3(:,3).*GAM(:,3);
    
    A =[dots(R2xR3,NT) + dots(S2,GV)...
        dots(R3xR1,NT) + dots(S3,GV)...
        dots(R1xR2,NT) + dots(S1,GV)];
    A =A./(asq*ones(1,3));
    SA(:,k)=A';
end

%************************************************************


function [keepProjPnt,keepTri,keeplambda,keepmu]=pointOnGeometry(VER,ITRI,ver)


minDist = inf;
keepProjPnt=[];
keeplambda=0;
keepmu=0;
keepTri=0;
for i=1:size(ITRI,1)
    vers = VER(ITRI(i,:),:);
    [result, projPnt, lambda, mu] = projectVerOnFace(vers,ver);
    if result ~= 0
        if ( lambda >= 0 &&...
                mu >= 0 &&...
                (lambda + mu) <= 1 )
            verOnFace = vers(1,:) * ( 1 - lambda - mu ) +  ...
                vers(2,:) * lambda +...
                vers(3,:) * mu;
        elseif ( lambda < 0 )
            lambda = 0;
            [~,verOnFace, mu]=edgeDist( ver, vers(1,:), vers(3,:) );
        elseif ( mu < 0 )
            mu = 0;
            [~,verOnFace, lambda]=edgeDist( ver, vers(1,:), vers(2,:) );
        else
            [~,verOnFace, mu]=edgeDist( ver, vers(2,:), vers(3,:) );
            lambda = 1.0 - mu;
        end
        d = norm( ver - verOnFace );
        if minDist >= d
            minDist = d;
            keepProjPnt = verOnFace;
            keepTri = i;
            keeplambda = lambda;
            keepmu = mu;
        end
    end
end


%%
function  [verOnline, d, alpha] = edgeDist(p, l1, l2)

%     // plindist(p, l1, l2, cp, lab): distance beween point p and
%     // line segment l1-l2. In cp the point on l1-l2 closed to
%     // p is returned; in lab the fraction of the distance
%     // l1-cp in respect to l1-l2 is returned.

s = ( l2 - l1 );
sq = sum(s.^2);

if ( sq == 0 )
    verOnline   = l1;
    alpha       = 0;
else
    alpha = min( max( 0.0, dot((p - l1),s) / sq), 1.0);
    verOnline   = l1 + s * alpha;
end
d = norm( p - verOnline );

%%


function [result, projPnt, lambda, mu] = projectVerOnFace(vers,pnt)

%     // Project the vertex on a face
result = 0;
lambda = 0;
mu     = 0;

s1 = ( vers(2,:) - vers(1,:) );
s2 = ( vers(3,:) - vers(1,:) );

n = [det([[1 0 0];s1;s2]) ...
    det([[0 1 0];s1;s2]) ...
    det([[0 0 1];s1;s2])];

tmp = sum(n.^2);

if ( tmp ~= 0 )
    %         /* To scale n such that is in a similar range as s1 and s2 and
    %          * consequently diff is in the same range as lambda and mu */
    tmp = sqrt( sqrt( tmp ) );
    n = n / tmp;
    
    sav = ( pnt - vers(1,:));
    det0 = det( [s1; s2; n] );
    if ( det0 ~= 0.0 )
        lambda = det( [sav; s2;  n] ) / det0;
        mu     = det( [s1;  sav; n] ) / det0;
        dif    = det( [s1;  s2; sav] ) / det0;
        if ( det0 < 0 )
            result = 1;
        else
            result = -1;
        end
    else
        error('point on ver error 2');
    end
else
    error('point on ver error 1');
end
projPnt = vers(1,:) + (vers(2,:)  - vers(1,:) ) * lambda + (vers(3,:)  - vers(1,:) ) * mu;




