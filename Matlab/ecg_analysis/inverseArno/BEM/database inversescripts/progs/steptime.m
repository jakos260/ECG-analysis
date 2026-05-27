for st=1:step,
itim=itim+1;
year=year+1;
survey(itim,1)=year;
pmi(1:101,2)=pmi(1:101,2)-mdeath;
pfi(1:101,2)=pfi(1:101,2)-fdeath;
pmi(101,2)=pmi(101,2)+pmi(100,2);
pfi(101,2)=pfi(101,2)+pfi(100,2);
for i=100:-1:2,
pmi(i,2)=pmi(i-1,2);
pfi(i,2)=pfi(i-1,2);
end
extra=aanwas;
pmi(1,2)=aanwas*fracmale;
pfi(1,2)=aanwas*(1.-fracmale);

ptot=pfi+pmi;
born=vrt(1:101,2).*pfi(1:101,2);
aanwas=sum(born);
mdeath=stm(1:101,2).*pmi(1:101,2);
fdeath=stf(1:101,2).*pfi(1:101,2);
exitmale=sum(mdeath);
exitfemale=sum(fdeath);
survey(itim,2)=sum(ptot(1:101,2));
survey(itim,3)=aanwas-exitmale-exitfemale;
survey(itim,4)=sum(ptot(21:66,2))/survey(itim,2);
aanwas-exitmale-exitfemale;
end
fourpl
saveasci('survey.lst',survey)