% normaliz.m
function NO=normaliz(A)
fac=norm(A,'fro');
NO=A/fac;
