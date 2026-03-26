% quadrangles.m
% function [quadrang,edge]=(VER,ITRI)
% compute the (internal) angles at the four vertices of the 
% quadrangle formed by two connected triangles 
% ITRI specifies three vertices for each if the TWO triangles
% VER specifies the vertex coordinates
% see also tri_angles


% 050201 a van oosterom

function [quadrang,edge]=quadrangles(VER,ITRI)
    QUAD=ITRI;
    edge=findedge(VER,QUAD,QUAD(1,1),3); 
    R(1,:)=VER(edge(2),:)-VER(edge(1),:);
    R(2,:)=VER(edge(3),:)-VER(edge(2),:);
    R(3,:)=VER(edge(4),:)-VER(edge(3),:);
    R(4,:)=VER(edge(1),:)-VER(edge(4),:);
    lr=norm3d(R);
    for i=1:4,
        j=icyc(i-1,4);, a=R(i,:); b=-R(j,:);
        sharp=sign(a*b');
        la=lr(i); lb=lr(j);
        C(i,:)=cross(a,b);
        acrb=norm3d(cross(a,b));
        sinang=asin(max(min(acrb/(la*lb),1),-1));
        if sharp<0, sinang=pi-sinang; end
        quadrang(i)=sinang;
             
    end
  % treat possibility of angles > pi 
  scc=sum(sign(C*C'));
  [mi imi]=min(scc);
  if mi<0, quadrang(imi)=2*pi-quadrang(imi);end 