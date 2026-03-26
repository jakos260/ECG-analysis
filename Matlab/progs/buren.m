% buren.m
% function [bver,BTRI]=buren(ITRI,node)
% for node: node of a triangulated surface,
%     find: 
%               bver=  indexes of the directly neighbouring vertices 
%               BTRI(:,1)   = triangle indices carrying the node
%               BTRI(:,2:4) = vertex indices of the triangle(s) carrying node
% A. van Oosterom
% NB: input arguments reduced; vectorized;  2004/10/30;
% node no longer included in bver
% see also: voisin, findedge, findloop, buurtris

% modified select; 2011-4-18

function [bver, BTRI]=buren(ITRI,node)

  bver=[]; BTRI=[];
  ntri=size(ITRI,1);
  

% use section below only when required while tuning
% itris=unique(ITRI); % all nodes referenced in ITRI
%   if ismember(node,itris)==0,
%      node
%      'warning: ismember(node,itris)==0'
%      pause
%   end
  
  select=find(ITRI(:,1)==node*ones(ntri,1)|ITRI(:,2)==node*ones(ntri,1)|ITRI(:,3)==node*ones(ntri,1));
  
  if size(select,1)==0, return, end
  
  BTRI=ITRI(select,:);
  %bver=skipreplica(BTRI(:));
  bver=unique(BTRI(:));
  BTRI=[select BTRI];
  bver(bver==node)=[];


