% fourier_contours.m
% NEWC=fourier_contours(CONT,kmax,fig,recs,contract);
% smoothing (set kmax) and resampling of nnodes of closed contours using
% Fourier expansion of (x,y,z) as a function of distance along the contour,
% CONT data may be closed but:
% reconstruction at ncrec DISJUNCT nodes
% if non-zero fig, process is monitored in figure(fig)
% contract<1 "pulls the string"; contract(>1 widens it; default:
% contract=1;

% A. van Oosterom; 20130724

function NEWC=fourier_contours(CONT,kmax,fig,recs,contract)

if nargin<5, contract=1; end

nct=size(CONT,1);
if sum(CONT(1,:)==CONT(nct,:),2)==3,
    CONT(nct,:)=[];
    nct=nct-1;
end

if nargin<=3,
    ncrec=nct;
end

if nargin>2,
    if fig~=0,
        figure(fig)
        clf
        hold on
        plot3(CONT(:,1),CONT(:,2),CONT(:,3),'*')
        grid on
        
    end
end

% path length along the contour
s=[0; cumsum(norm3d(CONT(2:nct,:)-CONT(1:nct-1,:)))];
smax=s(nct);

phi=s/smax*2*pi;
% NB: non-equal phi increments; and so:
% linear least squares approach to find Fourier expansion coefficients

if kmax >= floor(nct/4),
    kmax=floor(nct/4);
end

COSSIN=[ones(nct-1,1) cos(phi(1:nct-1)*(1:kmax)) sin(phi(1:nct-1)*(1:kmax))];
AB=pinv(COSSIN)*CONT(1:nct-1,:);  %(AB: Fourier expansion coeff)

% next: compute a fourier-interpolated contour

nint=401; % sets precission
phiint=(0:nint-1)'/(nint-1)*2*pi;

COSSINT=[ones(nint,1) contract*cos(phiint*(1:kmax)) contract*sin(phiint*(1:kmax))];
CONTINT=COSSINT*AB;

% find new path length
snew=[0; cumsum(norm3d(CONTINT(2:nint,:)-CONTINT(1:nint-1,:)))];
smaxnew=snew(nint);

% sub_sample fourier based curve at ncrec equally spaced intervals
% along the contour; 

if length(recs)==1,
    ncrecs=recs;
    dels=smaxnew/ncrecs;
else,
    ncrecs=length(recs);
end

if length(recs)>=2,
    % recs(i) = fraction of contour length
    select=[1 round(recs(2:end)*nint)]
else
    k=1;
    select=1;
    for i=2:ncrecs,
        level=(i-1)*dels;
        k=find( snew(1:nint-1)<=level & snew(2:nint)>level);
        select=[select,k];
    end 
end

NEWC=CONTINT(select,:);
if nargin>2,
    plot3(CONTINT(:,1), CONTINT(:,2), CONTINT(:,3),'m')
    title([' reconstruction based on kmax= ' num2str(kmax)])
    plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'r*')
end






