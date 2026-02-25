% tri_cleaning.m
% 20140409
% function [VER,ITRI,list]=tri_cleaning(VERIN,ITRIIN)
% removes any stray vertices from VER; 
% relabels ITRI accordingly
% list documents the original vertex labels

function [VER,ITRI,list]=tri_cleaning(VER,ITRI)
  vers=ITRI(:);
  vers=unique(vers); % all vertices named in ITRI
  nverintri=size(vers,1);
  nverin=size(VER,1);
  if nverin==nverintri, return,end
  
  LIST=[(1:nverin)' zeros(nverin,1)];
  LIST(vers,2)=1;
  VERLIST=[LIST(LIST(:,2)==1,1) (1:nverintri)'];
  VER=VER(VERLIST(:,1),:);
  INDEXLIST=[(1:nverin)' zeros(nverin,1)];
  INDEXLIST(VERLIST(:,1),2)=INDEXLIST(VERLIST(:,2));

  ITRI=[INDEXLIST(ITRI(:,1),2) INDEXLIST(ITRI(:,2),2) INDEXLIST(ITRI(:,3),2)]; 
  list=VERLIST(:,1);