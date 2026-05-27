% function parms = quamin_pvd(x,y,pinit,funtype)
%
% x         = length signal
% y     	= signal
% pinit     = initial estimates
% funtype   = sort of analysis [1 till 13], see rfunc.m
%
% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

function parms = quamin_pvd(x,y,pinit,funtype)
% script of: paramest
% Marquardt approach to solving a non-linear parameter estimation problem:
% See Avo - linsyst.pdf, p. 43-44 [section 6.2] for mathematical explanation.
 
% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute yest for initial estimate

[yest,G]    = rfunc(x,pinit,1,funtype);     % initial estimate [Y_est] based on pinit and function in funtype
res         = y - yest;                 	% residue between Y_meas and Y_est
normresiter = norm(res);                    % norm(residue)
iter        = 0;                            % start value iteration 
lambda      = 0;                            % start value lambda
temp        = [iter lambda normresiter];
ik          = 1;
RESNOW(ik,:)= temp;                         % output?
testnorm    = normresiter;                  % residue between measured and estimated T-wave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop1       = 1;
parms       = pinit;                        % initial parameters
npar        = length(pinit);                % number of initial parameters
if ~exist('noneg','var') == 1, noneg = 0; end
lambda      = 0.001;
maxlambda   = 1.e+10;
lambda_incr = 2;
lambda_decr = 4;
maxrelres   = 1.e-5;

% start iteration; outer loop
while loop1 >= 1,
   iter = iter + 1;             % increase k
   if iter > 100, break; end    % no more than 100 iterations for fitting
   
   GTG      = G'*G;             % compute Gradient matrix times Transpose of Gradient matrix
   gtgnorm  = norm(GTG,'fro');  % ??
   gtres    = G'*res;           % compute Gradient matrix times residue
   break1   = 0;
 
   % compute new estimate; if norm of residual does not decrease: 
   % restrict the step by constraining the step size, i.e. by increasing lambda
   while break1 == 0,
       lambda = lambda_incr*lambda;       % increase in lambda decreases magnitude of iteration step.
       if lambda > maxlambda, break1 = 1; break; end % if lambda is getting too big, break!
       
       % compute new estimate based on regularization parameter lambda;
       MAT  = GTG + lambda^2*eye(npar); % calculate left side Marquardt equation [AvO - linsyst.pdf [eq. (116)]]
       delp = pinv(MAT)*gtres;          % determine delta_p: [solving AvO - linsyst.pdf [eq. (116)]]
       
       % crude implementaion of noneg constraint [testp cannot go below zero] 
       if noneg == 1, testp = max(parms+delp',0); else testp = parms+delp'; end
       
       % test new estimate
       [yest]   = rfunc(x,testp,0,funtype); % new estimate with new initial parameters [mode = 0, no G output]
       res      = y - yest;
       testnorm = norm(res);
       sprintf('%0.4g  ',[iter lambda testnorm normresiter]);   % print [iteration number, lambda value, norm residue (lambda change: G same), norm residue (different G)]
       if testnorm < normresiter, break; end;
   end
    
    relres = abs(testnorm-normresiter)/normresiter;
    if relres < maxrelres, break1 = 1; end                  % check if residue is smaller than 1.e-5;
    if break1 == 1, break; end
    parms       = testp;                                    % set new parameter values, that decreased the residue between measured and estimated
    normresiter = testnorm;                                 % set new best residue for in lambda-loop
    lambda      = lambda / lambda_decr;                     % make lambda smaller
    [ytemp,G]   = rfunc(x,parms,1,funtype);                 % get new Gradient values
end
