
% pdisksub.m
% 'subroutine' variant of phidisk
% version: 20020730
% a. van Oosterom
% 
function phi=pdisksub(a,h,rho,kmax)

% using analytical expansion in Legendre plynomials

    hsq=h^2;
      r=sqrt(rho^2+hsq);
      mu=0; if r > 0, mu=h/r; end
      if r > a, phi=0;
        else
          phi=1.;
          if mu < 0, phi=-1; end
      end
      ar=0;
      if r > a, ar=a/r; end
      arsq=ar*ar;
      ra=r/a;
      rasq=ra*ra;
      faca1=-1*ra;
      faca2=0.5*arsq;

      for k=1:kmax,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for computing legendre polynomial of order zero and degree k
% note: use for k=1:kmax
if k==1,
   pk=mu;
   pkmm=1;
   pkm=mu;
else,
pk=((2*k-1)*mu*pkm-(k-1)*pkmm)/k;   
pkmm=pkm;
pkm=pk;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (rem is standard matlab remainder function)
        if rem(k,2)==1,
          if r > a,
             phi=phi+faca2*pk;
             faca2=-faca2*(k+2)/(k+3)*arsq;
          else,
             if r~=0, phi=phi+faca1*pk; end
             faca1=-faca1*k/(k+1)*rasq;
          end
          if mu==0, phi=0; end
        end
      end % end k-loop
  phi=phi/2;

