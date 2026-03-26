clf
[nr,nc]=size(survey);
plot(survey(:,1),survey(:,2),'white')
axis([survey(1,1) survey(nr,1) 0 2.e+7]);
hold on
plot(survey(:,1),2.e+7*survey(:,4),'yellow')
tekst=sprintf('%s','1.0');
text(survey(nr,1)+1,2.e+7,tekst)
tekst=sprintf('%s','0.5');
text(survey(nr,1)+1,1.e+7,tekst)
