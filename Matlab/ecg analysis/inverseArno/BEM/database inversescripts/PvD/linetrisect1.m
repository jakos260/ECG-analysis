% linetris.m
% function TRIINTER=linetris(VER,ITRI,l1,l2)
% find all intersections of linesegment running from location l1 to l2 with
% all elements a triangulated surface
% 20020116
% TRIINTER(:,1): triangle indexes of intersection 
% TRIINTER(:,2)=1 :view from outside, else: -1
% TRIINTER(:,3): (lambda) =fraction along edge 1
% TRIINTER(:,4): (mu)     =fraction along edge 2
% TRIINTER(:,5): (alpha)  =fraction along the line from l1 to l2
% no entry is included if the line runs parallel to the triangle

function TRIINTER=linetris(VER,ITRI,l1,l2)

TRIINTER=[];
dim=size(ITRI);
ntri=dim(1);
EPS=1e-9;
a1  =VER(ITRI(:,2),:)-VER(ITRI(:,1),:);
a2  =VER(ITRI(:,3),:)-VER(ITRI(:,1),:);
a3  =ones(ntri,1)*(l1-l2);
det0=det3d(a1,a2,a3);
index = 1:ntri;
A=[index' det0];
A=A(abs(det0)>eps,:);

k   =A(:,1);
l   =length(k);

det0=det0(k);
sav=ones(l,1)*l1-VER(ITRI(k,1),:);
det1=det3d(sav,a2(k,:),a3(k,:));
det2=det3d(a1(k,:),sav,a3(k,:));
det3=det3d(a1(k,:),a2(k,:),sav);
C=[k -sign(det0) det1./det0 det2./det0 det3./det0];
EPS=100000*eps;
TRIINTER=C(C(:,3)>=EPS&C(:,3)<=1-EPS&C(:,4)>=EPS&C(:,4)<=1-EPS&C(:,3)+C(:,4)<=1-EPS,:);
%TRIINTER=C(C(:,3)>=0&C(:,3)<=1&C(:,4)>=0&C(:,4)<=1&C(:,3)+C(:,4)<=1,:);


