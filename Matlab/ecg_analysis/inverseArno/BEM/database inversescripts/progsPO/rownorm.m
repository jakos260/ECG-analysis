function B=rownorm(A,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    rownorm
% 
% Desc:    applies norm to rows of a matrix and return the resut as a
% vector
% 
% Inputs:  A: input matrix. [opt]. option passed to norm. default is 2
% 
% Outputs: B: norm of rows of A
% 
% Created: 07-Nov-2013 11:16:18
% 
% Author:  peteroosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('opt','var')
    opt=2;
end
B=arrayfun(@(idx) norm(A(idx,:)), 1:size(A,1));



