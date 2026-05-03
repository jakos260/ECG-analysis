% getphi.m
% 20060508
% script of electrogr
% compute 1: infinite medium potentials at observation points
%         2: the potentials at the boundary
%         3: the potentials at the observation points
% for the homogeneous case: just 1) is required.

ral=a*sin(phia(t));
rbl=b*sin(phib(t));
yal=a*cos(phia(t));
ybl=b*cos(phib(t));

PHIOBS_inf=zeros(nn,ntims);
PHIB_inf=zeros(ntheta,ntims);

for k=2:ntims-1,
    VERS(1:nver/2,1:2)=VER(1:nver/2,1:2)*ral(k)/ainit;
    VERS(1:nver/2,3)=yal(k);
    if yal(k)==a,  VERS(1:nver/2,3)=b-rbl(k); end
    if yal(k)==-a, VERS(1:nver/2,3)=-b+rbl(k); end
    
    VERS(nver/2+1:nver,1:2)=VER(nver/2+1:nver,1:2)*rbl(k)/binit;
    VERS(nver/2+1:nver,3)=ybl(k);
    if ybl(k)==b, VERS(nver/2+1:nver,3)=a+ral(k);end
    
    
    % infinite medium potentials at OBS
    
    for i=1:nobs, 
       obs=OBS(i,:);
       [SA,index]=dsa(VERS,ITRI,obs,.2);% SA solid angles; units:ster-radians
       PHIOBS_inf(i,k)=sum(sum(SA))/(4*pi);
    end
    
    if inho==1,
       % infinite medium potentials at Boundary
       
       for i=1:ntheta,
           obs=a*VERB(sel(i),:);
           [SA,index]=dsa(VERS,ITRI,obs,.2);% SA solid angles; units:sterradians
           PHIB_inf(i,k)=sum(sum(SA))/(4*pi);
        end
    end   
end
PHIOBS=PHIOBS_inf;

if inho==1,
   % potentials at interface (boundary) with radius a
   PHIB=2*sigs/(sigm+sigp)*T*PHIB_inf; % T matrix as computed in 'electrog'
   %                                   % T=inv(I-kappa*B) accounts for sigm
                                                 
   % compute transfer from secondary sources at boundary to observation points
   W=zeros(1,ntheta);
   for i=1:nobs,
      obs=OBS(i,:);
      row=rowforw(a*VERB,ITRIB,obs);
      locate(i)=sum(row);
      W(1,1)=row(1);
      W(1,ntheta)=row(sel(ntheta));
      for jj=2:ntheta-1,
          W(1,jj)=sum(row(sel(jj):sel(jj+1)-1));
      end  
   
      % create the total solution at obs
   
      sigobs(i)=sigm; if locate(i) <.1, sigobs(i)=sigp; end
      if locate(i)>1.5 | locate(i) < .5,
         PHIOBS_sec(i,:)=  0.5*(sigm-sigp)/sigobs(i)*W*PHIB;
         PHIOBS(i,:)=sigs/sigobs(i)*PHIOBS_inf(i,:)+ PHIOBS_sec(i,:);
      else,% obs on boundary; assign value at nearest vertex
         dist=VERB(sel,:)-ones(ntheta,1)*[abs(obs(1)) obs(2:3)];
         [mi im]=min(norm3d(dist)');
         PHIOBS(i,:)=PHIB(im,:);
      end 
   end
end


