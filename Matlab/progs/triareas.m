% triareas.m
% function [verareas,triareas]=triareas(VER,ITRI)
% verareas: the areas around nodes VER
% triareas: the areas of individual triangles
% to compute total surface area, use: sum(...)
% see also: trinormals
% 050523 a van oosterom


function [verareas,triareas]=triareas(VER,ITRI)
   [ntri ndum]=size(ITRI);
   [nver ndum]=size(VER);
   triareas=zeros(ntri,1);
   verareas=zeros(nver,1);
for i=1:ntri,
    i1=ITRI(i,1);
    i2=ITRI(i,2);
    i3=ITRI(i,3);
    rm=VER(i2,:)-VER(i1,:);
    rp=VER(i3,:)-VER(i1,:);
    triareas(i)=norm(cross(rp,rm))/2;
    verareas([i1 i2 i3])=verareas([i1 i2 i3])+triareas(i)/3;
end







