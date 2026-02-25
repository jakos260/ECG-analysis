function [S,f,SNR]=noise_check(x,Fs,Fc)

[S,f]=pwelch(x,128,127,[],Fs);
h=abs(f-Fc);
[mn,i]=min(h);
SNR=20*log10(sum(S(1:i))/sum(S(i+1:end)));
