% tsvd.m
% function PHI=tsvd(PSI,k);   % uses singular values 1:k;  k=min([k ni,nj]);
% function PHI=tsvd(PSI,k,l); % uses singular values k>=1  ....  l <= min(sizePSI);

% A. van Oosterom, 20130425

function PHI=tsvd(PSI,k,l)
[ni,nj]=size(PSI);
[U,S,V]=svd(PSI);

if nargin==2,
    k=min([k ni,nj]);
    PHI=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
else,
    k=max(k,1);
    l=min([l,ni,nj]);
    k=min(k,l);
    PHI=U(:,k:l)*S(k:l,k:l)*V(:,k:l)';
end
     
