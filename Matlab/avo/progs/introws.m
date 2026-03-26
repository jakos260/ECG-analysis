function INTR=introws(M)
% integrate rows of matrix M
% cumsum type of action BUT: use local quadratic approx 
% [m n]=size(M); optimal if n is odd
% 20041016; A. van Oosterom

[m n]=size(M);
if n<2, 'number of samples should be greater than 1',return, end
INTR=zeros(m,n);

w=[5/12 8/12 -1/12];

for j=2:n-1;
   INTR(:,j)=INTR(:,j-1)+M(:,j-1:j+1)*w';
end

if rem(n,2)==0, 
    INTR(:,n)=INTR(:,n-1)+[M(:,n) M(:,n-1) M(:,n-2)]*w';
else,
    INTR(:,n)=INTR(:,n-1)+(M(:,n-1)+M(:,n))/2;
end
