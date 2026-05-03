% decompose_em
% 2006-08-30
% emperical mode decomposition (EMD)of signal
% decompose a signal by sucsessivley identifyig a signal passing through the mean of the
% two signals passing the maxima and the minima of the iterate, and subtracting it from
% the iterate
% function SIGCOMP=decompose_em(signal,ncomp,rdthresh,maxit)
% ncomp  the desired number of components
% rdthreshold relative difference used in terminating the identification of
% the individual components; default: 0.1
% maxnit restricting the maximumnumber of iterations when identifying the
% individual components default: 40

function SIGCOMP=decompose_em(signal,ncomp,rdthresh,maxit)
if nargin<4, maxit=40; end
if nargin<3, rdthresh=0.1; maxnit=40; end

nt=length(signal);
psi=signal; psi_prev=zeros(1,nt);
xx=1:nt;

SIGCOMP=[];
for icomp = 1:ncomp,
    
  for it=1:maxit,
      tmin=[1];
      tmax=[1];
      % find local extremes
      for i=2:nt-1,
          if psi(i)>psi(i-1) & psi(i)>psi(i+1),
             tmax=[tmax i];
          end   
          if psi(i)<psi(i-1)& psi(i)<psi(i+1),
             tmin=[tmin i];
          end
      end
      tmax=[tmax nt]; tmin=[tmin nt];
      ups=spline(tmax,psi(tmax),xx);
      downs=spline(tmin,psi(tmin),xx);
      comp=(ups+downs)/2;
      show=0;
      if show==1,
         figure(1)
         clf
         plot(psi,'k')
         hold on
         plot(ups,'r')
         %plot(tmax,psi(tmax),'r+')
         plot(downs,'b')
         %plot(tmin,psi(tmin),'b+')
         figure(2)
         clf
         plot(psi,'k') 
         hold on
         plot(psi-comp,'m')
         %pause
      end
      psi=psi-comp;
      rd=norm(comp)/norm(psi);
      ['reldiff with previous estimate of component ' num2str(icomp) ' :' num2str(rd)]
      
      if rd<rdthresh, break, end
  end
   
   SIGCOMP=[SIGCOMP; psi];
   psi = signal - sum(SIGCOMP,1);
   rchange=norm(psi-psi_prev)/norm(psi);
   ['relative change of including this component=' num2str(rchange)]
   psi_prev=psi;
end
   
