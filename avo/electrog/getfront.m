%getfront.m
% script of  electrog.m
% 20030128
% computes activation front
% on thick sphere myocardial model
% activation according to shortest myocardial path (Huygens)

if stim=='epi',
   % epicardial stimulus
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % first step: timing boundary points 
   % epicardial stimulus site
   for i=1:nphi
     if phi(i)<=2*phi1,
        dbf(i)=b*sqrt(2-2*cos(phi(i))); 
     else 
        dbf(i)=2*sqrt(b^2-a^2)+(phi(i)-2*phi1)*a;
     end
     if phi(i)<=phi1,
        daf(i)=sqrt(a^2+b^2-2*a*b*cos(phi(i))); 
     else 
        daf(i)=sqrt(b^2-a^2)+(phi(i)-phi1)*a;
     end
   end
   dbmax=max(dbf);
   dbmin=min(dbf);
   damax=max(daf);
   damin=min(daf);
   
   % find inverse values: angles for constant increments of the
   % distance functions
   phia=zeros(nphi,1);
   phib=zeros(nphi,1);
   for i=1:ntims,
      dist=(i-1)*dbmax/(ntims-1);
      
      % find intersection with dbf (local linear approx)
      for j=1:nphi-1
         if dbf(j)<=dist & dbf(j+1)>=dist,
            phib(i)=phi(j)+(dist-dbf(j))/(dbf(j+1)-dbf(j))*pi/(nphi-1);
         end
      end    

      % find intersection with abf (local linear approx)
      for j=1:nphi-1
          if daf(j)<=dist & daf(j+1)>=dist,
          phia(i)=phi(j)+(dist-daf(j))/(daf(j+1)-daf(j))*pi/(nphi-1);
          end
      end    
   end       

end

if stim=='endo',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endocardial stimulus

% first step: compute shortest path of boundary  points to 
% endocardial stimulus site
clf

% critical angle for un-obscured observation from stimulus site

for i=1:nphi

  if phi(i)<=phi1,
     dbf(i)=sqrt(b^2+a^2-2*a*b*cos(phi(i))); 
  else
     dbf(i)=a*(phi(i)-phi1)+sqrt(b^2-a^2);
  end
     daf(i)=a*(phi(i)); 
end

dbmax=max(dbf);
dbmin=min(dbf);
damax=max(daf);
damin=min(daf);

% find inverse values: angles for constant increments of the
% distance functions

phia=zeros(nphi,1);
phib=zeros(nphi,1);
for i=1:ntims,
  dist=(i-1)*dbmax/(ntims-1);

% find intersection with dbf (local linear approx)
  for j=1:nphi-1
      if dbf(j)<=dist && dbf(j+1)>=dist
        phib(i)=phi(j)+(dist-dbf(j))/(dbf(j+1)-dbf(j))*pi/(nphi-1);
        plot([0 phib(i) phib(i)],[dist dist 0],'k'), break,
      end
  end    

% find intersection with abf (local linear approx)
  for j=1:nphi-1
      if daf(j)<=dist && daf(j+1)>=dist,
      phia(i)=phi(j)+(dist-daf(j))/(daf(j+1)-daf(j))*pi/(nphi-1);
      plot([0 phia(i) phia(i)],[dist dist 0],'m'), break,
      end
  end    
end       

