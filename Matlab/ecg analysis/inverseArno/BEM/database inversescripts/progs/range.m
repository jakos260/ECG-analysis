% file range.m
% RANGE=range(X);
% computes the range of each column of X

% 20051201; A. van Oosterom

function RANGE=range(X);
ni=size(X);
if ni==1,
   RANGE=[min(X'); max(X')];
else
     RANGE=[min(X); max(X)];
 end



