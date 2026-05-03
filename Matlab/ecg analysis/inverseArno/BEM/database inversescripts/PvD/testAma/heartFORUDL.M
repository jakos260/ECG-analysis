	% forudl.m
% solving the neuman problem in matlab
% for a closed UDL (EDL) source
% Its geometry should be entered as the final entry in the geometry specs
clear
next =1;
% geometry specification


if next==1,
	ns=1;
	surf(1,:)=['uprofile.tri       '];  %torso
	% this final geometry should carry the UDL
	sigm(1)=1;
	outputname=['uprofile_avo.mat'];
end




nver=0;
ntri=0;
VER=[];
for is=1:ns,
   [VERS, ITRIS]=loadtri(surf(is,:));
%   if is==ns,refnode=VERS(251,:); end
   [nvers ndum]=size(VERS);
   [ntris ndum]=size(ITRIS);
   nver=nver+nvers;
   PNTSPE(is,:)=[nvers nver-nvers+1 nver];
   SIGMAS(is,2)=sigm(is);
   VER=[VER;VERS];
end

PNTSPE

% compute B matrix
B=zeros(nver,nver);
for js=1:ns,
	js
	[VERS,ITRIS]=loadtri(surf(js,:));
	for is=1:ns,
		for i=PNTSPE(is,2):PNTSPE(is,3),
			[B(i,PNTSPE(js,2):PNTSPE(js,3)),jsing]=rowforw(VERS,ITRIS,VER(i,:));
			k=jsing+PNTSPE(js,2)-1;
			if jsing~=0&js~=is, 
				B(i,k)=B(i,k)-1;
			end
		end
	end    
end

% compute surface containements
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

G=B(:,PNTSPE(ns,2):PNTSPE(ns,3));
nsrc=PNTSPE(ns,1);

for j=1:nsrc,
    k=nver-nsrc+j;
    for i=1:nver,
        if i==k, G(i,j)=G(i,j)-1; end
    end 
end

G=G*20;% NB: = 40 (mV) * G / 2  ( rowb=1/2pi SA)
[ndum nsrc]=size(G);

% scale by conductivities
for is=1:ns,
    for i=PNTSPE(is,2):PNTSPE(is,3),
        G(i,:)=2*G(i,:)/(SIGMAS(is,1)+SIGMAS(is,2));
    end
end
clear B
PSI=inv(eye(nver)-C)*G;
clear C
% reflate the PSI matrix

next=1;
if next==1,
for k=1:nsrc,
	sums=zeros(ns,1);
	for is=1:ns,
		if SIGMAS(is,1)~=0,
			sums(is)=sums(is)+sum(PSI(PNTSPE(is,2):PNTSPE(is,3),k));
		end  
	end
	for is=1:ns,
		theta=0;
		for js=1:ns,
			if SIGMAS(js,1)~=0,
				theta=theta+(SIGMAS(js,2)-...
					SIGMAS(js,1))*GAMMA(is,js)*sums(js)/(2*SIGMAS(js,1)*PNTSPE(js,1));
			end 
		end
		if SIGMAS(is,1)~=0,
			for i=PNTSPE(is,2):PNTSPE(is,3),
				PSI(i,k)=PSI(i,k)+theta;
			end
		end
	end   	  
end
end

% refer to WCT
%ref=mean([PSI(1,:);PSI(2,:);PSI(66,:)]);
%PSI=PSI-ones(nver,1)*ref;
size(PSI)
savemat(outputname,PSI);

