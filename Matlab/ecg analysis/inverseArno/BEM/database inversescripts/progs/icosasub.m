% icosasub.m
% generates icosahedron and, if desired, its refinements
% function [VER,ITRI]=icosasub(nref,mode)
% VER, ITRI : vertices and triangle indexes
% nref:       the number of times the icosaheder is to be refinded
% if mode==1, (default) vertex and triangle indices are ordered in clockwise fashion
%             as seen from the outside around node 1 (the north pole)

% matlab version: 2005-04-08; A. van Oosterom

function [VER,ITRI]=icosasub(nref,mode)

% nodes: basic isocahedron
ICO=[1 1 1 1 1 4  4  4  5  5  6  6 2 2  3  9 10 11  7   8 ;
     5 6 2 3 4 9 10  5 11  6  7  2 8 3  9 10 11  7  8   9 ; 
     4 5 6 2 3 3  9 10 10 11 11  7 7 8  8 12 12 12 12  12 ;];
ITRI=ICO';

VER( 1,:)=[0 0  1];
VER(12,:)=[0 0 -1];
rho=0.4*sqrt(5.);
fpi=.4*pi;
for i=2:11,
    im=i-2;
    if i>6, im=im-0.5; end
	VER(i,:)=[rho*cos(im*fpi) rho*sin(im*fpi) 0.5*rho];
     if i>6, VER(i,3)=-VER(i,3); end
end
    
if nref~=0,
   for refine=1:nref,
       ntri=size(ITRI,1);
       nver=size(VER ,1);
       TRI=[];
       for i=1:ntri,
           for j=1:3,
               jp=icyc(j+1,3);
               i1=ITRI(i,j); i2=ITRI(i,jp);
               h(1:3)=(VER(i1,:)+VER(i2,:))/2;
               h=h/norm(h);
               n=[];
               n=find(sum(VER==ones(nver,1)*h,2)==3);
               if isempty(n)==1,
                  VER=[VER;h];
                  nver=nver+1;
                  index(j)=nver;
               else,
                  index(j)=n;
               end
            end
            TRI=[TRI; 
            ITRI(i,1)  index(1)  index(3);
            index(1) ITRI(i,2)  index(2);
            index(3)  index(1)  index(2);
            index(3)  index(2) ITRI(i,3);];
        end
        ITRI=TRI;
   end
end

   if nargin>1&mode~=1, return,end
      nver=size(VER,1); ntri=size(ITRI,1);
 % relabel vertex indices; clockwise as seen from the outside, around node 1
      LIST=sortrows([(edgedist(nver,ITRI,1))' + (atan2(VER(:,1),VER(:,2))+pi)/(2*pi) (1:nver)'],1);
      list=[LIST(:,2) (1:nver)'];
      VER=VER(list(:,1),:);
      list=sortrows(list,1);
      ntri=size(ITRI,1);
      for i=1:ntri,
          ITRI(i,1:3)=[list(ITRI(i,1),2) list(ITRI(i,2),2) list(ITRI(i,3),2) ];
      end
      dist=(edgedist(nver,ITRI,1))';
      D=[dist(ITRI(:,1)) dist(ITRI(:,2)) dist(ITRI(:,3))];
      list=ceil(mean(D,2));
      CGRAV=(VER(ITRI(:,1),:)+VER(ITRI(:,2),:)+VER(ITRI(:,3),:))/3;
      LIST=sortrows([list + (atan2(CGRAV(:,1),CGRAV(:,2))+pi)/(2*pi) (1:ntri)'],1);
      %list=[LIST(:,2) (1:nver)'];
      %sort(LIST)
      ITRI=ITRI(LIST(:,2),:);
      

   