% gettdom
% function [tdom,ttm]=gettdom(PHI,jpoint)
% 2011-06-25
% Estimate dominant T wave
% from a weighted mean of T waves contained in PHI
% jpoint is a time marker after which depolarization is
% assumed to be completed
% ttm is timing apex Tdom
% exponential curve used in estimate during phase coinciding with QRS



function [tdom,ttm]=gettdompvd(PHI,qrsOnset,jpoint)

[nlds nt]=size(PHI);
Jampl = PHI(:,jpoint);
% PHI= lowpassma(PHI,50);
% qrsOnset = 200;
t= 1:jpoint-qrsOnset;
slope = -abs(Jampl) * 50 ./(jpoint - qrsOnset); 
% slope = Jampl ./(jpoint - qrsOnset); 
TQRS = t - (jpoint - qrsOnset)/2.0 ;
PSIQRS= (Jampl*ones(1,size(TQRS,2)))./(1+exp(slope*TQRS));	

% PSIQRS= slope * [0:(jpoint - qrsOnset)]; 
PSI = zeros(nlds,nt);
PSI(:,qrsOnset:jpoint-1) = PSIQRS;
PSI(:,jpoint:nt)=PHI(:,jpoint:nt);
PSI= lowpassma(PSI,40);

% estimate of tdom from weighted leads
tdom=ones(1,nt) * (PSI' * PSI);
ttm=find(tdom==max(tdom));

% final part adapted 20110623
% 
% next=0;
% if next==1,
% % find extrapolation down to t=1
% 
% tau=round((ttm-jpoint)*.3);
% t1=jpoint+20;
% y1=tdom(t1);
% 
% while y1<0,
%     t1=t1+1;
%     y1=tdom(t1);
% end
% 
% t2=t1+tau;
% t3=t2+tau;
% y1=tdom(t1);
% y2=tdom(t2);
% 
% if y2 < y1, y2=1.01*y1; end
% y3=tdom(t3);
% 
% alpha=(y3-y2)/(y2-y1);
% fac=alpha^(t1/tau);
% b=(y2-y1)/((alpha-1)*fac);
% a=y1-(y2-y1)^2/(y3-2*y2+y1);
% if alpha > 1,
%    if a < 0, 
%    a=0;
%    alpha=y2/y1;
%    b=y1/alpha^(t1/tau);
% end
% for j=1:t1-1,
%    tdom(j)=a+b*alpha^(j/tau);
% end 
% else,
%     
% end
    
    
% extrapolation down to t=1  

deriv=diffrows(tdom(jpoint:jpoint+100));
deriv=diffrows(tdom(qrsOnset:jpoint+100));
p=polyfit(0:length(deriv)-1,deriv,2);
yapprox = polyval(p,0:length(deriv)-1);
%  figure(1000);clf
% 
%  plot(yapprox)
%  hold on
%  plot(deriv,'r')
 a=find(yapprox == min(yapprox(1:round(length(yapprox)*0.7))));
 y1 = tdom(min(length(tdom),jpoint+a));
 y1pr = mean(yapprox(max(1,a-2):a+2));
 disp(['slope at jpoint ' num2str(y1pr)])
 
%  y1=tdom(jpoint);
%  y1pr=deriv(1);
 
%  % parabolic
%  ab=[2*jpoint   1;
%      jpoint^2    jpoint;]\[y1pr) y1]'
%     
%    for j=1:jpoint-1,
%        tdom(j)=ab(1)*j^2+ab(2)*j;
%    end 

% exponential
tau=y1/y1pr;
% y0=y1*exp(-jpoint/tau);
y0=y1*exp(-(jpoint+a)/tau);

j=1:jpoint+a;
tdom(j)=y0*exp(j/tau);
     
tdom=tdom/sum(tdom);