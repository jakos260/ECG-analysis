% file plotjes.m
function plotjes(PSI,labs,kleur)
dim=size(labs);
ni=dim(1);
nj=dim(2);
for i=1:ni,
for j=1:nj,
nplot=(i-1)*nj+j;
subplot(nj,ni,nplot)
k=labs(nplot);
plot(PSI(k,:),kleur);
tekst=sprintf(' lead: %d ',k);
title(tekst)
hold on
end
end
