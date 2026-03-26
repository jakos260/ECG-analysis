% linetris.m
% function TRIINTER=linetris(VER,ITRI,l1,l2)
% find all intersections of line directed from locations (vectors)l1 to l2 in 3D space
% with all triangles of a triangulated surface
% TRIINTER(:,1)   : triangle indices of all intersections 
% TRIINTER(:,2)=1 : line from l1 to l2 enters the triangle from the 
%                   outside domain, else: -1
% TRIINTER(:,3)   : (lambda) =fraction along edge 1
% TRIINTER(:,4)   : (mu)     =fraction along edge 2
% TRIINTER(:,5)   : (alpha)  =fraction along the line from l1 to l2
% TRIINTER=[] if the line runs parallel to the triangle

% NB: line 30: tolerance setting

% A. van Oosterom; 20140825


function TRIINTER=linetris(VER,ITRI,l1,l2)

TRIINTER=[];
ntri=size(ITRI,1);

a1  =VER(ITRI(:,2),:)-VER(ITRI(:,1),:);
a2  =VER(ITRI(:,3),:)-VER(ITRI(:,1),:);
a3  =ones(ntri,1)*(l1-l2);
det0=det3d(a1,a2,a3);
ind = 1:ntri;
A=[ind' det0];

k=A(abs(A(:,2))>1000*eps,1); % excludes triangles parallel to the line from l1 to l2

nk=size(k,1);

% find intersections of the planes containing the triangles selected and
% the line from l1 to l2 using Cramer's rule

det0=det0(k);
sav=ones(nk,1)*l1-VER(ITRI(k,1),:);

det1=det3d(sav    , a2(k,:), a3(k,:));
det2=det3d(a1(k,:), sav    , a3(k,:));
det3=det3d(a1(k,:), a2(k,:), sav    );

C=[k -sign(det0) det1./det0 det2./det0 det3./det0];

TRIINTER=C( (C(:,3)>=0 & C(:,3)<=1 &  ...
             C(:,4)>=0 & C(:,4)<=1 &  ...
             C(:,3)+ C(:,4)<=1 )   |  ...
             C(:,5)==0             |  ...
             C(:,5)==1,:);


