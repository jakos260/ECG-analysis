% fourier_contours.m
% NEWC=fourier_contours(CONT,kmax,fig,nrecs,contract);
% smoothing (set kmax) and resampling of nnodes of closed contour CONT(:,[x y z]) using
% Fourier expansions of CONT(:,[x,y,z)]) as a function of distance along the contour,
% CONT data need not be closed;
% reconstruction at ncrecs DISJUNCT nodes
% if non-zero fig, process is monitored in figure(fig)
% contract<1 "pulls the together "; contract(>1 widens it; default:
% contract=1;

% A. van Oosterom; 20131107

function NEWC=fourier_contours(CONT,kmax,fig,nrecs,contract)

if nargin<5, contract=1; end

nct=size(CONT,1);
if sum(CONT(1,:)==CONT(nct,:),2)==3,
    CONT(nct,:)=[];
    nct=nct-1;
end

if nargin<=3, 
    ncrecs=nct;
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

% resample
NEWC=zeros(nrecs,3);
NEWC(1,:)=CONTINT(1,:);

dels=smaxnew/nrecs;
for j=2:nrecs,
    t=(j-1)*dels;
    k=find(snew<=t);
    k=max(k);
    f=(t-snew(k))/(snew(k+1)-snew(k));
    NEWC(j,:)=(1-f)*CONTINT(k,:)+f*CONTINT(k+1,:); 
end

if nargin>2,
    if fig~0,
    plot3(CONTINT(:,1), CONTINT(:,2), CONTINT(:,3),'m')
    title([' reconstruction based on kmax= ' num2str(kmax)])
    plot3(NEWC(2:end,1),NEWC(2:end,2),NEWC(2:end,3),'k*')
    plot3(NEWC(1,1),NEWC(1,2),NEWC(1,3),'*r')
    end
end






