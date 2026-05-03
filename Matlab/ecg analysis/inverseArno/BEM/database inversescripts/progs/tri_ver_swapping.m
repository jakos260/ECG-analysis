% tri_ver_swapping.m
% 20131014
% [VER,ITRI]=tri_relabel(VER,ITRI,);
% vertices (RELABS(:,1) are swapped with those at locations: RELABS(:,2)
% ITRI indices are adapted accordingly

function [VER,ITRI]=tri_ver_swapping(VER,ITRI,SWAP);

% Test configuration
% [VER,ITRI]=make_sphere;
% figure(1)
% clf
% triplot
% x=lottery(12,12);
%  x=[10    11     2    12     7     1     4     6     9     8     3     5]'; 
% SWAP= [ (1:12)' x]  

nswaps=size(SWAP,1);

VERSAV=VER;

for i=1:nswaps,
    VER(SWAP(i,2),:)= VERSAV(SWAP(i,1),:);
end

ITRISAV=ITRI;

for i=1:nswaps,
    ITRI(ITRISAV==SWAP(i,1))=SWAP(i,2);
end

ITRSAV=[];

% figure(2)
% clf
% triplot

