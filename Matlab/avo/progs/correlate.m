function RHO=correlate(A)
dim=size(A);
RHO=zeros(dim(1),dim(2));
for i=1:dim(1),
for j=i+1:dim(2),
rho=corrcoef(A(i,:),A(j,:));
RHO(i,j)=rho(2,1);
end
end
