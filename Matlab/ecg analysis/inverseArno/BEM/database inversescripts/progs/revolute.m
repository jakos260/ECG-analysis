% function LIST=revolute(LIST,k)
% revolute the row elements of LIST k places;
% 0<=k<=size(LIST,2); 
% example:
% for LIST= [-3 -2 -1 0 1 2 3 4 5 6 7 ] , and k=2 it produces
%     LIST= [-1 0 1 2 3 4 5 6 7 -3 -2 ] 
% AvO; 201129

function LIST=revolute(LIST,k)
n=size(LIST,2);
LIST=LIST(:,([k+1:n 1:k]));