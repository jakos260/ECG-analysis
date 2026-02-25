% file contingency.m
% compute chi-square value of a contingenty table C

C=floor(C);
sumj=sum(C,2)
sumi=sum(C,1)
sumsum=sum(sumj);
CT=[C sumj;
    sumi sumsum]
    [ni nj]=size(C);
    df=(ni-1)*(nj-1)
chisq=0;
for i=1:ni,
    for j=1:nj,
        e=sumi(j)*sumj(i)/sumsum;
        chisq=chisq+(C(i,j)-e)^2/e;
    end
end
chisq


