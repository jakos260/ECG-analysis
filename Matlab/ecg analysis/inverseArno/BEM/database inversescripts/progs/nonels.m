% nonels.m
% NOn NEgative constrained linear Least Squares solution having minimal
% norm (NONEGLLSQMN)
% function [x,iter,res,ng]=nonels(E,f,xinit)

function [x,iter,res,ng]=nonels(E,f,xinit)

% 061104  A. van Oosterom; 
% solves the least squares problem Ex=f, with  x >= 0. , such that ||x|| is minimum
% compare: NNLS subroutine of Lawson and Hanson :
% 'Solving Least Square Problems' that also treats x>=0
% but does not return the minimal norm solution in the case of
% a rank deficient under-determined system of eqns; neither does the
% function lsqnonneg of MATLAB

% ng is the negative gradient vector; points downhill; available on output
% res is the vector of resiuals; available on output
% iter=number of iterations performed
% initial x: xinit can be specified, e.g. least squares minimum norm solution:
%       x=pinv(E)*f, x=max(0,x); default: x=0;
%       ia   indicates whether constraints are active:
%       ia()=1 : active constraint (at boundary, obviously!)
%       ia()=0 : constraint not active

% parameters setting the  
  tol=10*max(size(E))*norm(E,1)*eps;
  alpha=.1; maxalpha=1.e6;
  maxiter=300;
  log=0;  % log==1 will print various steps ofthe process 
  
 % initialize
  [m,n]=size(E);
  if nargin<3, xinit=zeros(n,1);, else, xinit=max(xinit,0); end
  
  prenor=inf;
  x=xinit;
  
  
% main loop
  brk=0;
 
  for iter=0:maxiter,
     % compute negative gradient vector: ng=-E'(Ex-f); points downhill
     % start by computing res=Ex-f
     res=E*x-f;
     resnorm=norm(res);
     if log==1,
     'iter/1000, x, resnorm, prenor, normx'
     [iter/1000 x' resnorm prenor norm(x)]
     end
     if resnorm>prenor;  brk=1; break; end
     prenor=resnorm;
     % now compute ng
     if iter >0,
        newg=-E'*res;
        %change=[];
        %change=ng(abs(ng-newg) >tol);
        ng=newg;
        %if isempty(change)==1, 'brk', brk=2, break, end
    else,
        ng=-E'*res;
     end
     % check Kuhn Tucker conditions
       ia=zeros(n,1);
       LIST=zeros(n,2); LIST(find(x<=-tol),1)=1; LIST(find(ng<tol),2)=1;
       ia(LIST(:,1)+LIST(:,2)>1)=1;
       % nact: number of active constraints
       nact=sum(ia);
       zg=[]; zg=ng(abs(ng)<=tol);
       % nzg:  number of zero gradients
       nzg=length(zg);
       [nact nzg];
       if log==1,
          'x, ng, LIST, ia:'
          [x ng LIST ia]
       end
       if nact+nzg >=n ; brk=3; break; end
     
     if brk>0; end
     % lsi: least squares with inequality constraints
         % copy the actively constrained columns of E to submatrix Ec,
         % the remainig ones to Eu,
         % the constrained elements of x to xc,
         % then find the LSMN solution to Eu*xu = f-Ec*xc
         jcon=find(ia==1);
         if isempty(jcon)==0,
             xc=x(jcon);
             Ec=E(:,jcon);
             fu=f-Ec*xc;
         else,
             fu=f;
         end
         
         juncon=find(ia==0);
         xu=x(juncon);
         Eu=E(:,juncon);
         restest1=norm(fu-Eu*xu);
         xutest=pinv(Eu)*fu;
         restest2=norm(fu-Eu*xutest);
         reldiftest=(restest1-restest2)/(restest1+restest2);
         if log==1,
            'xu, xutest ng(juncon)'
            [xu xutest ng(juncon)]
            'restest1, restest2 reldiftest'
            [restest1 restest2 reldiftest]
         end
         % xutest should be nonnegative; reldiftest should not be negative
         adapt=0;
         if (isempty(find(xutest<=-tol))==0) | reldiftest<0, adapt=1; end
         step=xutest-xu;
         nstep=norm(step);
         ngju=norm(ng(juncon));
         fac=norm(step)/norm(ng(juncon));
         ntel=0;
         while adapt==1,
            alpha=10*alpha;
            if log==1, alpha, end
            if alpha>maxalpha, brk=4; break, end
            ntel=ntel+1;
            if ntel>20, ntel, brk=5, break, end
            xutest=xu+(.01*step+alpha*ng(juncon)*fac)/(1+alpha^2);
            restest2=norm(fu-Eu*xutest);
            reldiftest=(restest1-restest2)/(restest1+restest2);
            if log==1,
               'xu, xutest ng(juncon)'
               [xu xutest ng(juncon)]
               'restest1, restest2 reldiftest'
               [restest1 restest2 reldiftest]
            %pause
         end
         
            if (reldiftest>=0)&isempty(find(xutest<=-tol))==1, adapt=0; end
         end
       
    % weave solution xutest into solution x; 
         alpha=sqrt(alpha);
         x(juncon)=xutest;
  end % end main loop
