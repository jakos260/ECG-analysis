% rfunc.m
% function [Yest,G]=rfunc(X,p,mode,functype);
% p=parameter vector
% mode==0: function values only; else function values and gradients
% functype==1   rflin0: straight line; zero intercept; requires 1 pinit value
%         ==2   rfpoly: general k-th order polynomial; requires k+1 pinit values 
%         ==3   rflogist: logist; 3 pinit values             
%         ==4   rfexpo: exponential+ shift; 3 pinit values 
%         ==5   rfexp0: exponential;        2 pinit values             
%         ==6   rflogistprod: product of two logistic functions exponential; 5 pinit values 
%         ==7   rflogistprod_g: product of two logistic functions exponential+Gauss; 8 pinit values 
%         ==8   rf_2expos:
%         ==9   rf_tangent: (not in progs)

%         ==10  rflogistprods:         
%         ==11  rflogist: logist + shift; 4 parms
%         ==12  rfgauss_over_sqrt; 4 parms
%         ==13  rflogist_prods



function [Yest,G]=rfunc(X,p,mode,functype);

 if functype==1;
    if mode==0, 
       [Yest]   =rflin0(X,p,mode);
    else
        [Yest,G]=rflin0(X,p,mode);
    end
 end
 
 if functype==2;
    if mode==0, 
       [Yest]  =rfpoly(X,p,mode);
    else
       [Yest,G]=rfpoly(X,p,mode);
    end
 end
 
if functype==3,
    if mode==0, 
       [Yest]  =rflogist(X,p,mode);
    else
       [Yest,G]=rflogist(X,p,mode);
    end
end
  
if functype==4,
   if mode==0,
      [Yest]  =rfexpo(X,p,mode);
   else
      [Yest,G]=rfexpo(X,p,mode);
   end
 end
  
if functype==5,
   if mode==0,
      [Yest]  =rfexp0(X,p,mode);
   else
      [Yest,G]=rfexp0(X,p,mode);
   end
end

if functype==6,
    if mode==0, 
       [Yest]  =rflogistprod(X,p,mode);
    else
       [Yest,G]=rflogistprod(X,p,mode);
    end
end

if functype==7,
    if mode==0, 
       [Yest]  =rflogistprod_g(X,p,mode);
    else
       [Yest,G]=rflogistprod_g(X,p,mode);
    end
end
    
if functype==8,
    if mode==0, 
       [Yest]  =rf_2expos(X,p,mode);
    else
       [Yest,G]=rf_2expos(X,p,mode);
    end
end

if functype==9,
    if mode==0, 
       [Yest]  =rf_tangent(X,p,mode);
    else
       [Yest,G]=rf_tangent(X,p,mode);
    end
end

if functype==10,
    if mode==0, 
       [Yest]  =rflogistprods(X,p,mode);
    else
       [Yest,G]=rflogistprods(X,p,mode);
    end
end

if functype==11,
    if mode==0, 
       [Yest]  =rflogist(X,p,mode);
    else
       [Yest,G]=rflogist(X,p,mode);
    end
end

if functype==12,
    if mode==0, 
       [Yest]  =rfgauss_over_sqrt(X,p,mode);
    else
       [Yest,G]=rfgauss_over_sqrt(X,p,mode);
    end
end

if functype==13,
    if mode==0, 
       [Yest]  =rflogist_prods(X,p,mode);
    else
       [Yest,G]=rflogist_prods(X,p,mode);
    end
end

