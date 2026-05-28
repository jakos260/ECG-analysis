% linetris.m
% function TRIINTER=linetris(VER,ITRI,l1,l2)
% find all intersections of linesegment running from locations (vectors)
% l1 to l2 in 3D space with all triangles of a triangulated surface
% TRIINTER(:,1)   : triangle indices of intersections 
% TRIINTER(:,2)=1 : line from l1 to l2 enters the triangle from the 
%                   outside domain, else: -1
% TRIINTER(:,3)   : (lambda) =fraction along edge 1
% TRIINTER(:,4)   : (mu)     =fraction along edge 2
% TRIINTER(:,5)   : (alpha)  =fraction along the line from l1 to l2
% TRIINTER=[] if the line runs parallel to the triangle
% A. van Oosterom; 20071012


function TRIINTER=linegeomsect(VER,ITRI,l1,l2,limit)

EPS=1e-11;
ntri=size(ITRI,1);

a1  =VER(ITRI(:,2),:)-VER(ITRI(:,1),:);
a2  =VER(ITRI(:,3),:)-VER(ITRI(:,1),:);
a3  =ones(ntri,1)*(l1-l2);
det0=det3d(a1,a2,a3);
ind = 1:ntri;
A=[ind' det0];

A=A(abs(A(:,2))>EPS,:);

k   =A(:,1);
l   =length(k);

det0=det0(k);
sav=ones(l,1)*l1-VER(ITRI(k,1),:);
det1=det3d(sav,a2(k,:),a3(k,:));
det2=det3d(a1(k,:),sav,a3(k,:));
det3=det3d(a1(k,:),a2(k,:),sav);
C=[k -sign(det0) det1./det0 det2./det0 det3./det0];
if limit
    TRIINTER=C(C(:,3) >= 0 & C(:,3) <= 1 & ...
               C(:,4) >= 0 & C(:,4) <= 1 & ...
               C(:,3) + C(:,4) <= 1 & ...
               C(:,5) > EPS & C(:,5) < 1-EPS,:);
    
else
    TRIINTER=C(C(:,3) >= 0 & C(:,3) <= 1 & ...
               C(:,4) >= 0 & C(:,4) <= 1 & ...
               C(:,3)+ C(:,4)<=1 & ...
               abs(C(:,5)) > EPS & abs(1-C(:,5)) > EPS,:);

end