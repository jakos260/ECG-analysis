function f=hazfun(par,tvals)
f= (1./(1.+par(1)/100*(par(2)-tvals))).^par(3);

