function i=icyc(k,n)
%function i=icyc(k,n)
i=k;
if k>n,
i=rem(k,n);
end
if k<=0,
i=rem(k,n)+n;
end
if i==0,
i=i+n;
end