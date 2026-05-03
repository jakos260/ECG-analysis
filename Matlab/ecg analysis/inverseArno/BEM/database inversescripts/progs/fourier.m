% [a0, a, b, kmax]=fourier(sig)
% find fourier coefficiens a0, a(k),b(k); k=1,kmax=(nt-1)/2 of signal: sig
% NOTE: fourier series expansion implies periodic extrapolation beyond the
% sampling interval
% frequencies are muliples of 1/(nt*sampling interval)=b1/(duration of
% the signal)

function [a0,a,b,kh]=fourier(sig)

% resynthesize for kh harmonics as:
  % sig(1:nt)=a0*ones(1,nt);
  % omega=2*pi/nt;
  % t=1:nt;
  % k=1:kh;
  % ANG=k'*omega*(t-1);
  % C(1:kh,1:nt)=cos(ANG);
  % S(1:kh,1:nt)=sin(ANG);
  % sig=sig+a*C+b*S;
  

dim=size(sig);
nt=max(dim);

F=fft(sig);
a0=real(F(1))/nt;
kh=floor((nt-1)/2);
a(1:kh)=2*real(F(2:kh+1))/nt;
b(1:kh)=-2*imag(F(2:kh+1))/nt;

