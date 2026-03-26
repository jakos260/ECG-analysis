% quamin.m
% script of: paramest
% Marquardt approach to solving a
% non-linear parameter estimation problem
% function values to be fitted should be present as column vector y
% absis values should be specified in calling script as column cector x
% initial parameters should be in rowvector pinit
% output=estimated parameters, contained in: parms

% A. van Ooserom; 2013_06_23 % notation changes implemented

% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute yest for initial estimate pinit
[yest,G]=rfunc(x,pinit,1,funtype);
funtype
res=y-yest;
normresiter=norm(res);
iter=0; lambda=0; ik=0;
temp=[iter lambda normresiter];
%'iter lambda normresiter'
%sprintf('%0.4g  ',temp)
ik=1;
RESNOW(ik,:)=temp;
testnorm=normresiter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outerloop=1;
parms=pinit;
npar=length(pinit);
if ~exist('noneg')==1, noneg=0;,end
lambda=0.001;

% start iteration; outer loop
while outerloop==1,
    %[iter normresiter]
    iter=iter+1;
    if iter>1000, break, end
    
    % compute GTG and G'*res
    GTG=G'*G;
    gtgnorm=norm(GTG,'fro');
    gtres =G'*res;
    
    % compute new estimate; if norm of residual does not decrease:
    % restrict the step by constraining the step size, i.e. by increasing
    % lambda
    innerloop=1;
    while innerloop==1,
        lambda=2*lambda;
        if lambda > 1.e+8, innerloop=0; lambda, bestnorm=testnorm; break, end
        % compute new estimate based on regularization parameter lambda;
        
%         lambda
%         GTG
%         npar
%         iter
        
        MAT=GTG+lambda^2*eye(npar);
       
        delp=pinv(MAT)*gtres;
        if noneg==1,
            % crude implementation of noneg constraint
            testp=max(parms+delp',0);
        else,
            testp=parms+delp';
        end
        
        % test new estimate
        [yest]=rfunc(x,testp,0,funtype);
        res=y-yest;
        testnorm=norm(res);
        sprintf('%0.4g  ',[iter lambda testnorm normresiter]);
        if testnorm < normresiter, break, end;
    end % end innerloop
    
    relres=abs(testnorm-normresiter)/normresiter;
%     if relres < 1.e-5,
%         sprintf('niter   relres  resnorm\n %3d   %0.4g   %0.4g',[iter relres testnorm]),
%     end
    if innerloop==0, break, end
    parms=testp;
    normresiter=testnorm;
    lambda=lambda/4;
    [yest,G]=rfunc(x,parms,1,funtype);
end % end outerloop

temp=[iter lambda normresiter];
'iter lambda normresiter'
sprintf('%0.4g  ',temp)

