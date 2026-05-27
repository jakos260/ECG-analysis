% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands
% rfunc.m
% function [Yest,G] = rfunc(X,p,mode,functype);
% p=parameter vector
% mode==0: function values only; else function values and gradients
% functype==1 : straight line; zero intercept; requires 1 pinit value
%         ==2   general k-th order polynomial; requires k+1 pinit values 
%         ==3   logist; requires 3 pinit values             
%         ==4   exponential+ shift; requires 3 pinit values 
%         ==5   exponential; requires 2 pinit values             
%         ==6   product of two logistic functions exponential; 3 pinit values 



function [Yest,G]=rfunc_pvd(X,p,mode,functype)

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
		Yest=rflogistprod_notch(X,p,mode);
	else
		[Yest,G]=rflogistprod_notch(X,p,mode);
    end
end



if functype==10,
    if mode==0, 
		Yest=rflogistprod_Tdom(X,p,mode);
	else
		[Yest,G]=rflogistprod_Tdom(X,p,mode);
    end
end

if functype==11,
    if mode==0, 
		Yest=rflogistprod_spike(X,p,mode);
	else
		[Yest,G]=rflogistprod_spike(X,p,mode);
    end
end


if functype==12,
    if mode==0, 
		Yest=rflogistprod_notch_only(X,p,mode);
	else
		[Yest,G]=rflogistprod_notch_only(X,p,mode);
    end
end

if functype==13,
    if mode==0, 
		Yest=rflogistprod_Tdom_noplat(X,p,mode);
	else
		[Yest,G]=rflogistprod_Tdom_noplat(X,p,mode);
    end
end


if functype==15,
    if mode==0, 
		Yest=rflogistprod_Tdom(X,p,mode);
	else
		[Yest,G]=rflogistprod_Tdom(X,p,mode);
    end
end

