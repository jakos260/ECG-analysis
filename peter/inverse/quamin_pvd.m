% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

function parms=quamin_pvd(x,y,pinit,funtype)
% script of: paramest
% Marquardt approach to solving a
% non-linear parameter estimation problem
 
  
% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute yest for initial estimate

[yest,G] = rfunc(x,pinit,1,funtype);

res = y - yest;
normresiter=norm(res);
iter=0; lambda=0;
temp=[iter lambda normresiter];
% disp('iter lambda normresiter')
% disp(sprintf('%0.4g  ',temp))
ik=1;
RESNOW(ik,:)=temp;
testnorm=normresiter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop1=1;
parms=pinit;
npar=length(pinit);
if ~exist('noneg','var')==1, 
    noneg=0;
end
lambda=.001;
% start iteration; outer loop
while loop1 >= 1
    %[iter normresiter]
   iter = iter+1;
   if iter > 100
       break;
   end
   % compute GTG and G'*res
   GTG = G'*G;
   gtgnorm = norm(GTG,'fro');
   gtres = G'*res;
   break1 = 0;
 
   % compute new estimate; if norm of residual does not decrease: 
   % restrict the step by constraining the step size, i.e. by increasing
   % lambda
   while break1 == 0
       lambda=2*lambda;
       if lambda > 1.e+10
           break1 = 1; 
           break;
       end
       % compute new estimate based on regularization parameter lambda;
       MAT = GTG + lambda^2*eye(npar);
%        if any(isinf(MAT)) | any(isnan(MAT))
%            sto=1;
%        end
       delp = pinv(MAT) * gtres;
       if noneg==1
          % crude implementaion of noneg constraint
          testp=max(parms+delp',0);
	   else
          testp=parms+delp';
       end
       % test new estimate
       [yest]=rfunc(x,testp,0,funtype);
       res = y - yest;
       testnorm=norm(res);
       sprintf('%0.4g  ',[iter lambda testnorm normresiter]);
       if testnorm < normresiter
           break
       end;
   end
    
    relres=abs(testnorm-normresiter)/normresiter;
    if relres < 1.e-5,
%          disp(sprintf('niter   relres\n %3d   %0.4g',[iter relres]))
         break1=1;
    end
    if break1==1
        break
    end
    parms=testp;
    normresiter=testnorm;
    lambda = lambda / 4;
    [~,G] = rfunc(x,parms,1,funtype);
end




















