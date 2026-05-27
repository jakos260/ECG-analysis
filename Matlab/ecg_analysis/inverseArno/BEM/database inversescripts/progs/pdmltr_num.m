function [ phi ] = pdmltr_num( VER, obs )

%   Based on "TwoD" by Lawrence F. Shampine.
%   Ref: L.F. Shampine, "Matlab Program for Quadrature in 2D",
%   Appl. Math. Comp., 202 (2008) 266-274.


r1 = VER(1,:);
r2 = VER(2,:);
r3 = VER(3,:);
r = obs;

n = -cross(r2-r1,r3-r1);
A = norm(n);
%n = n/A;
phi(1) = A * quad2d(@fun1,0,1,0,@(x)(1-x));
phi(2) = A * quad2d(@fun2,0,1,0,@(x)(1-x));
phi(3) = A * quad2d(@fun3,0,1,0,@(x)(1-x));
phi = phi(:);

% -------------------------------------------------------------------------

function z = fun1(x,y)
    [row,col] = size(x);
    R = ones(numel(x),1)*(r1-r) + x(:)*(r2-r1) + y(:)*(r3-r1);
    Rnorm = sqrt(sum(R.^2,2));
    z = reshape(1./Rnorm .* [1-x(:)-y(:)],row,col);
end

function z = fun2(x,y)
    [row,col] = size(x);
    R = ones(numel(x),1)*(r1-r) + x(:)*(r2-r1) + y(:)*(r3-r1);
    Rnorm = sqrt(sum(R.^2,2));
    z = reshape(1./Rnorm .* x(:),row,col);
end

function z = fun3(x,y)
    [row,col] = size(x);
    R = ones(numel(x),1)*(r1-r) + x(:)*(r2-r1) + y(:)*(r3-r1);
    Rnorm = sqrt(sum(R.^2,2));
    z = reshape(1./Rnorm .* y(:),row,col);
end


end
