% skipreplicas.m
% SINGLES=skipreplicas(LIST,column,small)
% skipreplicas skips all rows in which column is replicated
% default column: 1
% default skipping criterium: max(abs(entries))<small
% A. van Oosterom; 2012/02/15
% note: the MATLAB script unique performs sorting, which is frequently not
% desired

function SINGLES=skipreplicas(LIST,column,small)
  if nargin==1, column=1; end
  if nargin<=2,small=eps; end 
 
  SINGLES=[];
  
  while ~isempty(LIST),
      k= abs(LIST(:,column)-LIST(1,column))<=small;
      SINGLES=[SINGLES;LIST(1,:)];
      LIST(k,:)=[];
  end
  
     
    
  
  
  
  
