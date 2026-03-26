% tri_relabel.m
% 20141015
% [VER,ITRI]=tri_relabel(VERIN,ITRIN,oldlabs,newlabs);
% labs pertain to VER rows; rows of VER are adapted accordingly
function [VER,ITRI]=tri_relabel(VERIN,ITRIN,oldlabs, newlabs);


noldlabs=max(size(unique(oldlabs)));
nnewlabs=max(size(unique(newlabs)));

if noldlabs~=nnewlabs, 'lab lists should be of equal length', pause, return, end

VER=zeros(size(VERIN,1),3);
for i=1:noldlabs,
    VER(newlabs(i),:)=VERIN(oldlabs(i),:);
end    
ITRIN;
ntris=size(ITRIN,1);
ITRI=zeros(ntris,3);

for i=1:ntris,
    for j=1:3,
        ITRI(i,j)=find(oldlabs==ITRIN(i,j));
    end
end       