% extremes.m
% function extr=extremes(A);
% extr=[min(min(A)) row_min column_min max(max(A)) row_max col_max]
function extr=extremes(A);

colmins = min(A);
[colmin j]=min(colmins);
[mini i]  = min(A(:,j));

colmaxs=max(A);
[colmax l] = max(colmaxs);
[maxi k]  = max(A(:,l));
extr=[mini i j maxi k l];

