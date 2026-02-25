% tri_ver_swapping.m
% 20131120
% [VER,ITRI]=tri_relabel(VER,ITRI,SWAP);
% vertices (SWAP(:,1) are swapped with those at locations: SWAP(:,2)
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

VER(SWAP(:,1),:)= VERSAV(SWAP(:,2),:);
VER(SWAP(:,2),:)= VERSAV(SWAP(:,1),:);
 
ITRISAV=ITRI;

for i=1:nswaps,
    ITRI(ITRISAV==SWAP(i,1))=SWAP(i,2);
    % NBNB next line was found to be required in some applications
    % needs to be explained
    %ITRI(ITRISAV==SWAP(i,2))=SWAP(i,1);
end

ITRSAV=[];

% figure(2)
% clf
% triplot

