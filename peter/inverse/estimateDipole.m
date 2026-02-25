function dip = estimateDipole(GEOM,ECG,useleads,dipLoc)


PSI = forwardDipole(GEOM,dipLoc);
PHI_dips = PSI(1:length(GEOM.TVER),:);
PHI_dips = PHI_dips - ones(length(PHI_dips ),1)*mean(PHI_dips(1:3,:));

L = PHI_dips(useleads,:);

% L=zeromean(PHI_dips(4:12,:)'); % lead vectors
% L=L';


% direct pseudo inverse
D = L \ ECG;
dip = D;

% % damped least squares
% mu=1.e-2;
%D=inv(L'*L+mu^2*eye(size(L,2)))*L'*PHImid;

% 
% % truncated SVD inverse
% [U,S,V]=svd(L);
% 
% k=min(9,size(L,2)); % inverse based on based on  k singular vector pairs
% 
% s_inv=diag(S);
% s_inv=1./s_inv(1:k);
% 
% T=V(:,1:k)*diag(s_inv)*(U(:,1:k)');
% % T=T';
% % T=zeromean(T);
% % T=T';

%D=T*PHImid;

PHI2 = L*D;




function PSI = forwardDipole(GEOM,dipLoc)


% dipLoc = mean(GEOM.VER);
DIPS = [dipLoc 1 0 0;
        dipLoc 0 1 0;
        dipLoc 0 0 1;];
nrsc = size(DIPS,1);    

nvers=0;
ntri=0;
VER = GEOM.TVER;
ITRI= GEOM.TITRI;
ns = 1;
nver = length(GEOM.TVER);
PNTSPE(1,:)=[nvers nver-nvers+1 nver];

SIGMAS = [1];

% compute B matrix
B=zeros(nver,nver);
for js=1:ns,
%     [VERS,ITRIS]=loadtri(surf(js).name);
    VERS= GEOM.TVER;
    ITRIS = GEOM.TITRI;
    for is=1:ns,
        for i=PNTSPE(is,2):PNTSPE(is,3),
            [B(i,PNTSPE(js,2):PNTSPE(js,3)),jsing] = rowforw(VERS,ITRIS,VER(i,:));
            k=jsing + PNTSPE(js,2)-1;
            if jsing~=0 && js~=is, 
                B(i,k) = B(i,k) - 1; 
            end
        end
    end
end

% compute surface containments: SURCON(i,j)= 1 for i==j, else,
%                                          = 2 if Sj contains Si; else
%                                          = 0
SURCON=eye(ns);
for is=1:ns,
    i1=PNTSPE(is,2);
    for js=1:ns,
        if js~=is,
            i2=PNTSPE(js,2);
            i3=PNTSPE(js,3);
            test=sum(B(i1,i2:i3));
            if test > 1.5, SURCON(is,js)=2; end
        end
    end
end

% SURCON
%pause
% compute deflation factors
[SIGMAS,KAPPA,DEFL,GAMMA]=deflat(SURCON,SIGMAS);

% deflate the B matrix
C=zeros(nver,nver);
for is=1:ns,
    for i=PNTSPE(is,2):PNTSPE(is,3),
        for js=1:ns,
            for j=PNTSPE(js,2):PNTSPE(js,3),
                C(i,j)=B(i,j)*KAPPA(is,js)-DEFL(is,js)/PNTSPE(js,1);
            end
        end
    end
end


% compute the infinite medium potentials in a medium of unit conductivity
% from dipoles;  position and strength:
% GG=zeros(nver,nsrc);
GG = cgmadipole(VER,DIPS);

G = GG;

% scale by conductivities
for is=1:ns,
    for i=PNTSPE(is,2):PNTSPE(is,3),
        G(i,:) = 2*G(i,:)/(SIGMAS(is,1)+SIGMAS(is,2));
    end
end

%savemat(outputname,inv(eye(nver)-C));

%'done'
%pause

PSI = inv(eye(nver)-C) * G;

% BMA=inv(eye(nver)-C);
% savemat([prefix '_bma.mat'],BMA)

next=1;
if next==1, % reflate the PSI matrix
    
    for k=1:size(DIPS,1)
        sums=zeros(ns,1);
        for is=1:ns,
            if SIGMAS(is,2)~=0,
                sums(is)=sums(is)+sum(PSI(PNTSPE(is,2):PNTSPE(is,3),k));
            end
        end
        for is=1:ns,
            theta=0;
            for js=1:ns,
                if SIGMAS(js,2)~=0,
                    theta=theta+(SIGMAS(js,1)-...
                        SIGMAS(js,2))*GAMMA(is,js)*sums(js)/(2*SIGMAS(js,2)*PNTSPE(js,1));
                end
            end
            if SIGMAS(is,2)~=0,
                for i=PNTSPE(is,2):PNTSPE(is,3),
                    PSI(i,k)=PSI(i,k)+theta;
                end
            end
        end
    end
end

%%
function G=cgmadipole(OBS,DIPS)
[nobs idum]=size(OBS);
[nsrc idum]=size(DIPS);       
G=zeros(nobs,nsrc);
for k=1:nsrc,
   R=OBS-ones(nobs,1)*DIPS(k,1:3);
   r=norm3d(R);
   ss=find(r>eps);  
   if isempty(ss)==0, 
       G(ss,k)=dots(R(ss,:),DIPS(k,4:6))./r(ss).^3;
   end
end

G = G/(4*pi);
