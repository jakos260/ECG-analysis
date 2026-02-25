% dml.m
% function  PHI=dml(VER,ITRI,obs)
% vectorized  version of pdmltr_anal.m
% compute the matrix PHI (size: ntri,3) of the contributions of the three distributed monolayers over each of
% the ntri triangles ITRI (each with unit strength at one vertex and zero strength at the other two vertices) 
% to the potential at observation point obs. 
% =vectorized  version of pdmltr_anal.m

% calling: cross; norm3d; det3d; dots; solida;
% % if index==0: obs does NOT lie in the plane of the triangle
% % index=1: obs is an external point of the triangle
% % index=2: obs lies on the edge of the triangle;
% % index=3: obs is an internal point of the triangle
% % index=4: obs is a node of the triangle

% A. van Oosterom; 2012_02_17

% NB:NB;NB;NB does not include scaling by sigma and 4*pi


function  PHI=dml(VER,ITRI,obs)
nver=size(VER,1);
ntri=size(ITRI,1);

VER=VER-ones(nver,1)*obs;

R1=VER(ITRI(:,1),:); 
R2=VER(ITRI(:,2),:);
R3=VER(ITRI(:,3),:);
LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];

% cross products ; normals to the tetrahedron faces other than the source
% triangle
Z1=cross(R2,R3);   % normal to the face not including vertex 1, etc
Z2=cross(R3,R1);   %  
Z3=cross(R1,R2);   %  

E1=R2-R1;  % edges 1 of all triangles; size(ntri,3)
E2=R3-R2;  %       2
E3=R1-R3;  %       3

LE=[norm3d(E1) norm3d(E2) norm3d(E3)];  %edge lengths; size(ntri,3)

E1n=E1./(LE(:,1)*ones(1,3));
E2n=E2./(LE(:,2)*ones(1,3));
E3n=E3./(LE(:,3)*ones(1,3)); % normalized edges


% computes solid angles (see solida as used by pdmltr_anal)
% blockproduct=triple vector product
block=dots(R1,Z1);
DOTS=[dots(R1,R2) dots(R2,R3) dots(R3,R1)];
denom=LR(:,1).*LR(:,2).*LR(:,3) + ...
      LR(:,1).*DOTS(:,2) + LR(:,2).*DOTS(:,3) + LR(:,3).*DOTS(:,1);
                                           
% index(abs(block)<=eps & denom >eps)                 =1;
% index(abs(block)<=eps & denom >= -eps & denom <=eps)=2;
% index(abs(block)<=eps & denom <  -eps)              =3;
% index(LR(:,1)<=eps|LR(:,2)<=eps|LR(:,3)<=eps)       =4;
 
sa=-2*atan2(block,denom)';   % solid angles subtended at obs by the triangles
sa(abs(block)<=eps)=0;

N=cross(E1,E2);
normn=norm3d(N);          % twice the area of the source triangle
Nu=N./(normn*ones(1,3));  % unit normal of the source triangle
h_tetra=block./normn;      % signed distance between [obs] and the
                          % plane of the source triangle 
                          % height of the tetrahedron 

% rows of WW: the signed distances between the projection of obs
%              on the plane of the source triangle and the edges 

WW=[det3d(R1,R2,Nu)./LE(:,1)    det3d(R2,R3,Nu)./LE(:,2)   det3d(R3,R1,Nu)./LE(:,3) ]; % dimension m^2

GAM=zeros(ntri,3);
a=LR(:,1); b=LR(:,2); c=LE(:,1);
GAM(:,1)=log(((b+c).^2 - a.^2 +ones(ntri,1)*eps)./(abs( b.^2-(a-c).^2 +ones(ntri,1)*eps)));

a=LR(:,2); b=LR(:,3); c=LE(:,2);
GAM(:,2)=log(((b+c).^2 - a.^2 +ones(ntri,1)*eps)./(abs( b.^2-(a-c).^2 +ones(ntri,1)*eps)));

a=LR(:,3); b=LR(:,1); c=LE(:,3);
GAM(:,3)=log(((b+c).^2 - a.^2 +ones(ntri,1)*eps)./(abs( b.^2-(a-c).^2 +ones(ntri,1)*eps)));
% size(ntri,3) ;rows: line integrals over edges of the source triangles

uniform= h_tetra .* sa' +  sum(WW.*GAM,2); % uniform=int(1/r) overe the triangle (obs is the origin)

Q= [LR(:,2).^2-LR(:,1).^2-LE(:,1).^2    LR(:,3).^2-LR(:,2).^2-LE(:,2).^2      LR(:,1).^2-LR(:,3).^2-LE(:,3).^2 ];

D= [Q(:,1).^2-4*LR(:,1).^2.*LE(:,1).^2  Q(:,2).^2-4*LR(:,2).^2.*LE(:,2).^2    Q(:,3).^2-4*LR(:,3).^2.*LE(:,3).^2];

I= LR(:,[2 3 1]).*LE/2  +  Q.*(LR(:,[2 3 1])-LR)./(4*LE) -  D./(8*LE.^2).*GAM; 
        
    M1=[dot(E2,E1n,2) dot(E3,E1n,2) dot(E1,E1n,2)];
    M2=[dot(E2,E2n,2) dot(E3,E2n,2) dot(E1,E2n,2)];
    M3=[dot(E2,E3n,2) dot(E3,E3n,2) dot(E1,E3n,2)];
    
    TEST= uniform*ones(1,3).*[sum(Z1.*Nu,2) sum(Z2.*Nu,2) sum(Z3.*Nu,2)]; 
      % = Z*u_n'
    ADD=-[ dot([M1(:,1) M2(:,1) M3(:,1)],I,2) dot([M1(:,2) M2(:,2) M3(:,2)],I,2) dot([M1(:,3) M2(:,3) M3(:,3)],I,2)];
      %= -  E(l,:)*En(k,:)' *I
       
    PHI=( TEST+ ADD)./(normn*ones(1,3));
  
    
% above:
% vectorized version of pdmltr_anal.m, which treats a single
% triangle;
% below: relevant lines of pdmltr_anal.m, which treats a single triangle;
%    k=[1 2 3];
%    l=[2 3 1]; 
%    m=[3 1 2];
%    h_tetra=block/normn;
%    ww(k)  = det3d(R(k,:),R(l,:),ones(3,1)*u_n)./e(k); % dimension m^2
%    gam=log ((r(l).*e(k)  + sum(R(l,:).*E(k,:),2)+ ones(3,1)*eps )./ ...
%            (abs(r(k).*e(k)  + sum(R(k,:).*E(k,:),2))+ ones(3,1)*eps ) );
%    uniform= h_tetra * omega +  ww'*gam; %dimension [m]  uniform (monolayer!!           
%    Q_sq= r(l).^2-r(k).^2-e(k).^2
%    D= Q_sq.^2   -4*r(k).^2.*e(k).^2
%    I=contour integral along the edges of the function r(lambda)
%    I= r(l).*e(k)/2  +  Q_sq.*(r(l)-r(k))./(4*e(k)) -  D./(8*e(k).^2).*gam;
%    phi=(Z*u_n' *uniform -  E(l,:)*En(k,:)' *I)/normn;
    
   




   
   
  
