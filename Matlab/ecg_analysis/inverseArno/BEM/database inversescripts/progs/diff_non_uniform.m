function deriv=diff_non_uniform(y,x)
% function deriv=diff_non_uniform(y,x)
% computing derivative of y based on of samples 
% UN-equally distributed along x-axis
% local parabole based estimate; 
    [nix njx]=size(x);
    if nix < njx ; x=x'; [nix njx]=size(x); end
    [niy njy]=size(y);
    if niy<njy; y=y'; [niy njy]=size(y); end
    if niy~=nix, 'incorrect input', pause, end
    n=nix;
    deriv=zeros(n,1);
    x3=x(1:3);
    y3=y(1:3);
    M=[x3.^2 x3 ones(3,1)];
    abc=pinv(M)*y3;
    deriv(1)=2*abc(1)*x3(1)+abc(2);

for j=2:n-1,
    x3=x(j-1:j+1);
    y3=y(j-1:j+1);
    M=[x3.^2 x3 ones(3,1)];
    abc=pinv(M)*y3;
    deriv(j)=2*abc(1)*x3(2)+abc(2);
end
    x3=x(n-2:n);
    y3=y(n-2:n);
    M=[x3.^2 x3 ones(3,1)];
    abc=pinv(M)*y3;
    deriv(n)=2*abc(1)*x3(3)+abc(2);