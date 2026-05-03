% testfitnodes3d.m
% Marquardt approach to solving a
% non-linear parameter estimation problem applied to fitting two sets
% of points in space by means of rotations followed by a shift
clear; 
datum=4.0306;
clf

% case settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load input data

Yref=[0 0 0; 1 0 0; 0 1 0; 0 0 1];

VER=Yref;
VALS=VER;
ITRI=[1 4 2;2 4 3;1 4 3;1 3 2];
triplot
txt=text(VER(:,1),VER(:,2),VER(:,3),'*','fontsize',16);

pause
X=rotash(Yref,[45/180 45/180 0.2],[0 .2 0]);
%X=X*1.1
set(hs,'Vertices',X);
pause


VER=X;
np=size(VER,1);
gravityy=mean(Yref)
gravityx=mean(X)


%X=X+ones(np,1)*(gravityy-gravityx);
%[mean(X);mean(Yref)]
%set(hs,'Vertices',X);
%pause

y=Yref(:);
pinit=zeros(1,6);
npar=length(pinit);

noneg=0;

  
% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute yest for initial estimate
[yest,G]=rf3dpnts(X,pinit,1);
res=y-yest;
normresiter=norm(res);
iter=0; lambda=0; ik=0;
temp=[iter lambda normresiter];
'iter lambda normresiter'
sprintf('%0.4g  ',temp)
ik=1;
RESNOW(ik,:)=temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loop1=1;
parms=pinit;
GTG=G'*G;
lambda=.001;

% start iteration; outer loop

while loop1>=1,
   iter=iter+1;
   
   % compute GTG and G'*res
   GTG=G'*G;
   gtgnorm=norm(GTG,'fro');
   gtres =G'*res;
   break1=0;
 
   % compute new estimate; if norm of residual does not decrease: 
   % restrict the step by constrainig the step size, i. e. by increasing lambda
   % 
   while break1==0,
       lambda=2*lambda;
       if lambda > 1.e+10*gtgnorm, break1=1; lambda,   bestnorm=testnorm;break, end
       % compute new estimate based on regularization parameter lambda;
       MAT=GTG+lambda^2*eye(npar);
       delp=pinv(MAT)*gtres;
       if noneg==1,
          % crude implementaion of noneg constraint
          testp=max(parms+delp',0);
       else,
          testp=parms+delp';
       end
       % test new estimate
       [yest]=rf3dpnts(X,testp,0);
       res=y-yest;
       testnorm=norm(res);
       'iter lambda testnorm normresiter'
        sprintf('%0.4g  ',[iter lambda testnorm normresiter])
       % pause
       if testnorm < normresiter, break, end;
    end
    
    if abs(testnorm-normresiter)/normresiter< 1.e-4, break1=1; end
    
    if break1==1, break, end
    
    parms=testp
    
    normresiter=testnorm;
    lambda=lambda/4;
    [yest,G]=rf3dpnts(X,parms,1);
    VER=reshape(yest,np,3);
    set(hs,'Vertices',VER);
pause

end

Yest=reshape(yest,np,3)
Yref

