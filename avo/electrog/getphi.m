% getphi.m
% 20030208
% script of electrogr
% compute 1: infinite medium potentials at observation points
%         2: the potentials at the boundary
%         3: the potentials at the observation points
% for the homogeneous case: just 1) is required.

% infinite medium potntials at Boundary
PHIOBS=zeros(nn,ntims);

% infinite medium potntials at OBS
PHIBinf=zeros(ntheta,ntims);

ral=a*sin(phia(t));
rbl=b*sin(phib(t));
yal=a*cos(phia(t));
ybl=b*cos(phib(t));

for k=2:ntims-1,
    VERS(1:nver/2,1:2)=VER(1:nver/2,1:2)*ral(k)/ainit;
    VERS(1:nver/2,3)=yal(k);
    if yal(k)==a,  VERS(1:nver/2,3)=b-rbl(k); end
    if yal(k)==-a, VERS(1:nver/2,3)=-b+rbl(k); end
    
    VERS(nver/2+1:nver,1:2)=VER(nver/2+1:nver,1:2)*rbl(k)/binit;
    VERS(nver/2+1:nver,3)=ybl(k);
    if ybl(k)==b, VERS(nver/2+1:nver,3)=a+ral(k);end
    
    for i=1:nn,
       obs=OBS(i,:);
       [SA,index]=dsa(VERS,ITRI,obs,.2);% SA solid angles; units:sterradians
       PHIOBS(i,k)=sum(sum(SA))/(4*pi);
    end
    
    if inho==1,
       for i=1:ntheta,
           obs=a*VERB(sel(i),:);
           [SA,index]=dsa(VERS,ITRI,obs,.2);% SA solid angles; units:sterradians
           PHIBinf(i,k)=sum(sum(SA))/(4*pi);
        end
    end   
end
 
if inho==1,
   % boundary potentials
   PHIB=2*sigs/(sigm+sigp)*T*PHIBinf;
    
   % compute transfer from secondary sources at boundary to observation points
   W=zeros(1,ntheta);
   for i=1:nn,
      obs=OBS(i,:);
      row=rowforw(a*VERB,ITRIB,obs);
      locate=sum(row);
      W(1,1)=row(1);
      W(1,ntheta)=row(sel(ntheta));
      for jj=2:ntheta-1,
          W(1,jj)=sum(row(sel(jj):sel(jj+1)-1));
      end  
      
      % create the total solution at obs
      sigobs=sigm; if locate <.1, sigobs=sigp; end
      if locate>1.5 | locate < .5,
         PHIOBS(i,:)=sigs/sigobs*PHIOBS(i,:)+0.5*(sigm-sigp)/sigobs*W*PHIB;
      else,% obs on boundary; assign value at nearest vertex
         dist=VERB(sel,:)-ones(ntheta,1)*[abs(obs(1)) obs(2:3)];
         [mi im]=min(norm3d(dist)');
         PHIOBS(i,:)=PHIB(im,:);
      end 
   end
end


