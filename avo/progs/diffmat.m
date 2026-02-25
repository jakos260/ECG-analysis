function D=diffmat(n)
for i=1:n,
D(i,1)=0;
end
D(1,1)=-.5;
D(2,1)=.5;
for j=2:n-1,
for i=1:n,
D(i,j)=0.;
end
D(j-1,j)=-.5;
D(j+1,j)=.5;
end
for i=1:n,
D(i,n)=0;
end
D(n-1,n)=-.5;
D(n,n)=.5;

