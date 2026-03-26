% tsvdinv.m
% function X=tsvd(A,B,k)
% solves A X = B as truncated inverse 
% k is truncation rank

function X=tsvd(A,B,k)
[U,S,V]=svd(A);

Sinv=inv(S(1:k,1:k));

X=V(:,1:k)*Sinv*U(:,1:k)'*B;
