function varargout = computeECGsimEGM(varargin)

sigmaG =  1;		% general
sigmaL =  0.2;	% lungs
sigmaLi = 0.4;    % liver
sigmaB =  3;		% blood
sigmaF =  0.16;  % fat
sigmaR =  0.06;  % ribcage
sigmaH =  1;     % heart= source
lp=1;

% sigmaG = 1;		% general
% sigmaL = 1;	% lungs
% sigmaLi= 1;    % liver
% sigmaB = 1;		% blood
% sigmaF = 1;  % fat
% sigmaR = 0.1;  % ribcage
% sigmaH = 1;     % heart= source


if nargin<4
    error('Not enough parameters!!!');
else
    GEOM = varargin{1};
    VOL = varargin{2};
    obsP=varargin{3};
    doplot=0;
    S=varargin{4};
    SAH = S;
    % if length(varargin)>4
    %     doplot=varargin{5};
    % end
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
                case 'liver'
                    sigmaLi=varargin{pp+1};pp=pp+2;
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

[AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'atria',obsP); 
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'lcav',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'rcav',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'ventr',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'rlung',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'llung',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'thorax',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'ribcage',obsP); end

if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'fatpad_1',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'fatpad_2',obsP); end
if isempty(AMA) [AMA, cavityStr] = computeAMAOnGeometry(GEOM,VOL,'liver',obsP); end

if ~isempty(AMA)
    disp(['location: ' cavityStr]);
    PHI = AMA * S;
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
insideLI = 0;
if sum(solida(VOL.geom.VER,VOL.geom.ITRI,obsP)) > 0.01
    % Within the source is undefined
    PHI = zeros(1,size(S,2));
    varargout{1} = PHI;
    if nargout > 1
        varargout{2} = 'in myocardium';
    end
    return;
elseif sum(solida(GEOM.thorax.VER,GEOM.thorax.ITRI,obsP)) < 0.01
    % Outside the body
    PHI = zeros(1,size(S,2));
    varargout{1} = PHI;
    if nargout > 1
        varargout{2} = 'on thorax';
    end
    return;
elseif isfield(GEOM,'rlung') && sum(solida(GEOM.rlung.VER,GEOM.rlung.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaL;
    cavityStr='right lung';
    insideRL = 1;
elseif isfield(GEOM,'llung') && sum(solida(GEOM.llung.VER,GEOM.llung.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaL;
    cavityStr='left lung';
    insideLL = 1;
elseif isfield(GEOM,'rcav') && sum(solida(GEOM.rcav.VER,GEOM.rcav.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaB;
    cavityStr='right cavity';
    insideRC = 1;
elseif isfield(GEOM,'lcav') && sum(solida(GEOM.lcav.VER,GEOM.lcav.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaB;
    cavityStr='left cavity';
    insideLC = 1;
elseif isfield(GEOM,'liver') && sum(solida(GEOM.liver.VER,GEOM.liver.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaLi;
    cavityStr='liver';
    insideLI = 1;
elseif isfield(GEOM,'ribcage') && sum(solida(GEOM.ribcage.VER,GEOM.ribcage.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaR;
    cavityStr='ribcage';
    insideRB = 1;
elseif isfield(GEOM,'fatpad_1') && ...
        sum(solida(GEOM.fatpad_1.VER,GEOM.fatpad_1.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaF;
    cavityStr='fatpad 1';
    insideF1 = 1;
elseif isfield(GEOM,'fatpad_2') && ...
        sum(solida(GEOM.fatpad_2.VER,GEOM.fatpad_2.ITRI,obsP)) > 0.01
    sigmaObs  = sigmaF;
    cavityStr='fatpad 2';
    insideF2 = 1;
end
% disp(['location in : ' cavityStr]);
if doplot, figure(12345);clf; hold on; ncomp=1;end

PHI = (40 * sigmaH / sigmaObs) * sa(VOL.geom.VER,VOL.geom.ITRI,obsP) * S;

if doplot, plot(lowpassma(PHI ,lp),'b','linewidth',1);lt{ncomp}='heart';ncomp = ncomp+1;end

if isfield(VOL,'THORAX')
    PHIT = ((sigmaG-0) / sigmaObs) * sa(GEOM.thorax.VER,GEOM.thorax.ITRI,obsP) * (VOL.THORAX *SAH);
    if doplot, plot(lowpassma(PHIT ,lp),'g','linewidth',1);lt{ncomp}= 'torso';ncomp = ncomp+1;end
    PHI = PHI + PHIT; % by definition + because always inside
end

if isfield(VOL,'RCAV')
    PHIRCAV = ((sigmaB-sigmaG) / sigmaObs ) * (sa(GEOM.rcav.VER,GEOM.rcav.ITRI,obsP) * (VOL.RCAV * SAH));
    if doplot, plot(lowpassma(PHIRCAV ,lp),'r','linewidth',1);lt{ncomp}= 'Rcav';ncomp = ncomp+1;end
    PHI=PHI + PHIRCAV;
end
if isfield(VOL,'LCAV')
    PHILCAV = ((sigmaB-sigmaG) / sigmaObs ) * sa(GEOM.lcav.VER,GEOM.lcav.ITRI,obsP) * (VOL.LCAV*SAH);
    if doplot, plot(lowpassma(PHILCAV ,lp),'r:','linewidth',1);lt{ncomp}= 'Lcav';ncomp = ncomp+1;end
    PHI=PHI + PHILCAV;
end
if isfield(VOL,'LIVER')
    PHILIV = ((sigmaLi-sigmaG) / sigmaObs ) * sa(GEOM.liver.VER,GEOM.liver.ITRI,obsP) * (VOL.LIVER*SAH);
    if doplot, plot(lowpassma(PHILIV ,lp),'k:','linewidth',1);lt{ncomp}= 'liver';ncomp = ncomp+1;end
    PHI=PHI + PHILIV;
end

if isfield(VOL,'RLUNG')
    PHIRL= ((sigmaL-sigmaG) / sigmaObs )*sa(GEOM.rlung.VER,GEOM.rlung.ITRI,obsP) * (VOL.RLUNG *SAH);
    if doplot, plot(lowpassma(PHIRL ,lp),'m','linewidth',1);lt{ncomp}= 'Rlung';ncomp = ncomp+1;end
    PHI = PHI + PHIRL;
end

if isfield(VOL,'LLUNG')
    PHILL = ((sigmaL-sigmaG) / sigmaObs ) * sa(GEOM.llung.VER,GEOM.llung.ITRI,obsP)*(VOL.LLUNG*SAH); 
    if doplot, plot(lowpassma(PHILL,lp),'m:','linewidth',1);lt{ncomp}= 'Llung';ncomp = ncomp+1;end
    PHI = PHI + PHILL;
end

if isfield(VOL,'RIBCAGE')
    PHIRIB = ((sigmaR-sigmaG) / sigmaObs ) * sa(GEOM.ribcage.VER, GEOM.ribcage.ITRI, obsP ) * ( VOL.RIBCAGE * SAH );
    if doplot, plot(lowpassma(PHIRIB ,lp),'c','linewidth',1);lt{ncomp}= 'ribca';ncomp = ncomp+1;end
    PHI = PHI + PHIRIB;    
end
if isfield(VOL,'FATPAD_1')
    PHIFAT1 = ((sigmaF-sigmaG) / sigmaObs ) * sa(GEOM.fatpad_1.VER, GEOM.fatpad_1.ITRI, obsP )*(VOL.FATPAD_1*SAH);
    if doplot, plot(lowpassma(PHIFAT1 ,lp),'y-','linewidth',1);lt{ncomp}= 'fat_1';ncomp = ncomp+1;end
    PHI = PHI + PHIFAT1;
end
if isfield(VOL,'FATPAD_2')
    PHIFAT2 = ((sigmaF-sigmaG)  / sigmaObs ) * sa(GEOM.fatpad_2.VER, GEOM.fatpad_2.ITRI, obsP )*(VOL.FATPAD_2*SAH);   
    if doplot, plot(lowpassma(PHIFAT2,lp),'y:','linewidth',1);lt{ncomp}= 'fat_2';ncomp = ncomp+1;end
    PHI = PHI + PHIFAT2;
end

if ~isempty(S) && doplot
    plot(lowpassma(PHI ,lp),'k','linewidth',1); lt{ncomp}= ['final in ' cavityStr] ; legend(lt);drawnow; 
    disp(['location in : ' cavityStr '  '  num2str(PHI)]);
end
varargout{1} = PHI;
if nargout > 1
    varargout{2} = cavityStr;
end
%************************************************************
function [AMA,cavityStr] = computeAMAOnGeometry(GEOM,VOL,geom,obsP)
hasGeometry = isfield(GEOM, geom);
cavityStr='';

AMA=[];
if hasGeometry
    small = 1e-5;
    eval(['A=GEOM.' geom ';'])
if strcmp(geom,'ventr')
    B = VOL.VENTRICLES;
elseif isfield(VOL, upper(geom))
    eval(['B=VOL.' upper(geom) ';'])
else
    return;   % <-- se non esiste il campo, esci senza errore
end
    % eval(['A=GEOM.' geom ';'])
    % if ( strcmp(geom,'ventr'))
    %     eval(['B=VOL.VENTRICLES;'])
    % else
    %     eval(['B=VOL.' upper(geom) ';'])
    % end

    [vProjPnt,vTri,vLambda,vMu]=pointOnGeometry( A.VER, A.ITRI,obsP );
    if norm(vProjPnt-obsP) < small
        itri= A.ITRI(vTri,:);
        cavityStr=['on ' geom];
        AMA = B(itri(1),:) + (B(itri(2),:) -B(itri(1),:)) * vLambda +...
                (B(itri(3),:) -B(itri(1),:)) * vMu;

        % if ( strcmp(geom,'ventr'))
        %     AMA = VOL.VENTRICLES(itri(1),:) + (VOL.VENTRICLES(itri(2),:) -VOL.VENTRICLES(itri(1),:)) * vLambda +...
        %         (VOL.VENTRICLES(itri(3),:) -VOL.VENTRICLES(itri(1),:)) * vMu;
        % else
        % 
        % end
    end
end

%************************************************************

function rowb=sa(VER,ITRI,obsP)

      
% OMEGA = distributedSa(VER,ITRI,obsP,0.01);
% rowb = zeros(1,length(VER));
% for j=1:length(ITRI),
%     rowb(ITRI(j,:)) = rowb(ITRI(j,:)) + OMEGA(:,j)';
% end
% rowb=rowb/(4*pi);

[rowb,jsing] = rowforw(VER,ITRI,obsP);

if jsing~=0,
    rowb(jsing)= rowb(jsing)+1;
end
rowb = rowb/2;
%%==========================================================

function [rowb,jsing]=rowforw(VER,ITRI,obs)

nver=size(VER,1);
ntri=size(ITRI,1);

rowb=zeros(1,nver);
OMEGA=distributedSa(VER,ITRI,obs,.01); % NOTE: uses dsa

for j=1:ntri,
    ij=ITRI(j,:);
    rowb(ij)=rowb(ij) + OMEGA(:,j)';
end

% search singularity

jsing=find(norm3d(VER-ones(nver,1)*obs)/norm(VER,'fro')<1.e-6);

if isempty(jsing), jsing=0; end

if jsing~=0,% treat singularity (determines auto solid angle)
    
    loop=loopnode(ITRI,jsing);
    
    nb=size(loop,2);
    for kk=1:nb,
        k=icyc(kk-1,nb); l=kk; m=icyc(kk+1,nb);
        sa(kk)=solida([VER(loop(k),:);VER(loop(l),:);VER(loop(m),:)],[1 3 2],obs);
    end
    ndpos=sum(sa>0);
    ndneg=sum(sa<0);
    
    if ndpos==nb || ndneg==nb % use spherical cap approximation
        % but only if theta, the cone angle of the cap as viewed from the origin
        % of the approximating sphere, is less than pi/6.
        % center of gravity of direct neighbours
        center=mean(VER(loop,:));
        
        % find rho: the radius of circle through direct neighbours
        % use: rms value of distances of direct neighbours to center
        rho=norm(VER(loop,:)-ones(nb,1)*center,'fro')/sqrt(nb);
        
        % distance from center to node(jsing)
        node2c=norm(VER(jsing,:)-center);
        coshalft=rho/sqrt(rho^2+node2c^2);
        theta=2*acos(coshalft);
        
        if theta < 2*pi/3
            if theta < 1.e-4,
                rowb(jsing)=pi*theta/2;
            else
                rowb(jsing)=4*pi*(1-coshalft)/theta;
            end
            if ndneg==nb, rowb(jsing)=-rowb(jsing);end
            rowb(loop)=rowb(loop)+(2*pi-sum(rowb))/nb;
        else
            rowb(jsing)=2*pi-sum(rowb);
        end
    else
        rowb(jsing)=2*pi-sum(rowb);
    end
end

rowb=rowb/(2*pi);

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




