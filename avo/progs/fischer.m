% fischer.m
% Fisher's linear discrimination ; two groups (classes)
% specify: GR1 (size(n1,nfeatures)) and GR2 (size(n2:,nfeatures))
%          f1 and f2: strings for features 1 and 2 :xlabel and ylabel        
% AvO, 080224

nn=size(GR1);
n1=nn(1);
nn=size(GR2);
n2=nn(1);
CV1=cov(GR1);
CV2=cov(GR2);
m1=mean(GR1);
m2=mean(GR2);
SW=CV1+CV2;

w=inv(SW)*(m1-m2)';
w=w/norm(w);                              % unit vector along line of separation

% m1a=(w'*m1')*w;
% m2a=(w'*m2')*w;

vn(1)=-w(2);
vn(2)=w(1);

c=(m1*n1+m2*n2)/(n1+n2);                  % weighted mean of group averages; weights: no of obs
d=c*w;

dextr=d;
scoremax=[0 0 0];
% scan along 
for dd=d*0.98:d/250:d*1.02, 
    ncorrect1=0;
    for i=1:n1,   
        test=GR1(i,:)*w-dd;
        if test >= 0 
           ncorrect1=ncorrect1+1;
        end
    end
    ncorrect2=0;
    for i=1:n2,
        test=GR2(i,:)*w-dd;
        if test <= 0
           ncorrect2=ncorrect2+1;
        end
    end
    score=[ncorrect1/n1 ncorrect2/n2 (ncorrect1+ncorrect2)/(n1+n2)];
   if score(3) > scoremax(3)
      scoremax=score;
      dextr=dd;
      scor=score;
   end
end

% mark class1 in scatterplot of features 1 and 2
plot(GR1(:,1),GR1(:,2),'red+')
hold on
% mark class2 in scatterplot of features 1 and 2
plot(GR2(:,1),GR2(:,2),'blue+')
xlabel(f1) % feature 1
ylabel(f2) % feature 2

plot([m1(1) m2(1)],[m1(2) m2(2)],'green') % vector normal to separation line/plane
scale=axis;
xmi=scale(1);
xma=scale(1);
ymi=scale(3);
yma=scale(3);
xtest=(dextr-scale(3)*w(2))/w(1); 
if(xtest > scale(1))
   xmi=xtest;
   if(xtest > scale(2)) xmi=scale(2);  end
end
xtest=(dextr-scale(4)*w(2))/w(1); 
if(xtest > scale(1))
    xma=xtest;
    if(xtest > scale(2)) xma=scale(2); end
end

ytest=(dextr-xmi*w(1))/w(2); 
if(ytest > scale(3)) 
    ymi=ytest;
    if(ytest > scale(4)) ymi=scale(4); end
end
ytest=(dextr-xma*w(1))/w(2); 
if(ytest > scale(3))
    yma=ytest;
    if(ytest > scale(4)) yma=scale(4); end
end

plot([xmi xma],[ymi yma],'r')  % line/plane of separation

ab=scale(1);
ae=scale(2);
bb=scale(3);
be=scale(4);

xtekst=ab- 0.15*(ae-ab);
ytekst=be+ 0.05*(be-bb);

tekst=sprintf('%s',class1,' n= ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab);
tekst=sprintf('%0.5g',n1);
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%s',' correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab);
tekst=sprintf('%0.5g',scor(1));
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.2*(ae-ab);
tekst=sprintf('%s',class2,'   n= ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab);
tekst=sprintf('%0.5g',n2);
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%s',' correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab);
tekst=sprintf('%0.5g',scor(2));
text(xtekst,ytekst,tekst);

xtekst=xtekst-0.24*(ae-ab);
ytekst=bb-.1*(be-bb);
tekst=sprintf('%s','total correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.24*(ae-ab);
tekst=sprintf('%0.5g',scor(3));
text(xtekst,ytekst,tekst);

xtekst=ab-0.15*(ae-ab);
tekst=sprintf('%s','w;d');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%0.5g',w(1));
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab);
tekst=sprintf('%0.5g',w(2));
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.2*(ae-ab);
tekst=sprintf('%0.5g', d);
text(xtekst,ytekst,tekst);

score=scor;