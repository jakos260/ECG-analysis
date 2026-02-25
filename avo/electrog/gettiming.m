% gettiming.m
% script of electrog
% 20030128
% computes timing of the passage of the wave front

% epicardial stimulus; get tims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stim=='epi ',
   % first step: timing at boundary points 
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
end

if stim=='endo',
   % endocardial stimulus
   % first step: compute shortest path of boundary  points to 
   for i=1:nphi
      if phi(i)<=phi1,
         dbf(i)=sqrt(b^2+a^2-2*a*b*cos(phi(i))); 
      else,
        dbf(i)=a*(phi(i)-phi1)+sqrt(b^2-a^2);
      end
         daf(i)=a*(phi(i)); 
   end
end

dbmax=max(dbf);
dbmin=min(dbf);
damax=max(daf);
damin=min(daf);
maxdist=max([damax dbmax]);

% next, find inverse values: phi values for constant increments of the
% distance functions
phia=zeros(ntims,1);
phib=zeros(ntims,1);
for i=1:ntims,
   dist=(i-1)*pi/(ntims-1);

   % find intersection with dbf (local linear approx)
   if dist>=dbmax,phib(i)=pi;else
      for j=1:nphi-1
         if dbf(j)<=dist & dbf(j+1)>=dist,
         phib(i)=phi(j)+(dist-dbf(j))/(dbf(j+1)-dbf(j))*pi/(nphi-1);
         , break,
         end
      end
   end    

   % find intersection with daf (local linear approx)
   if dist>=damax, phia(i)=pi;else,
   for j=1:nphi-1
       if daf(j)<=dist & daf(j+1)>=dist,
       phia(i)=phi(j)+(dist-daf(j))/(daf(j+1)-daf(j))*pi/(nphi-1);
       , break,
       end
      end
   end    
end       



