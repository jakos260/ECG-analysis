function m=std1d(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    std1d
% 
% Desc:    calculate std of x, ignoring NaN's
% 
% Inputs:  x: n-dimensianl matrix, possibly containing NaN's
% 
% Outputs: m: std of all elements of x not NaN.
% 
% Created: 12-Nov-2013 15:08:04
% 
% Author:  peteroosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(:);
m=std(x(~isnan(x)));