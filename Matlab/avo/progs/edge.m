function edge(ITRI,VER,range)
% plot the edge of a 2D plot of (a part of) a triangulated surface
% the range of triangles to be included is specified by range; 
% default: [1:ntri]

dim=size(ITRI);
ntri=dim(1);

if nargin==2, range=[1 ntri]; end
zin=zeros(ntri,1);
for i=range(1):range(2),
    
i1=ITRI(i,1);
i2=ITRI(i,2);
i3=ITRI(i,3);
r1=VER(i2,:)-VER(i1,:);
r2=VER(i3,:)-VER(i1,:);
r3=cross(r1,r2);
if r3(1) > 0, zin(i)=1;  end;
if r3(1)== 0, zin(i)=0;  end;
if r3(1) < 0, zin(i)=-1; end;
end

% spot and plot edge 
for i=range(1):range(2)-1,
   tritsi=ITRI(i,:);
   loop=0;
   for j=i+1:range(2),
      if zin(i) ~= zin(j),
        tritsj=[ITRI(j,3) ITRI(j,2) ITRI(j,1)];
        k=0;
        while k<3,
           k=k+1;
           i1=tritsi(k);
           l=0;
           while l<3,
              l=l+1;
              j1=tritsj(l);
              if i1==j1, 
                  m=k+1; if m>3, m=1; end
                  i2=tritsi(m);
                  n=l+1; if n>3, n=1; end
                  j2=tritsj(n);
                  if i2==j2,
                      plot3([VER(i1,1) VER(i2,1)],[VER(i1,2) VER(i2,2)],[VER(i1,3) VER(i2,3)],'-k')
                      loop=1; break,
                  end
             end
             if loop==1, break, end
             end
             if loop==1, break, end
          end
          if loop==1, break, end
      end
      if loop==1, break, end
   end
end
