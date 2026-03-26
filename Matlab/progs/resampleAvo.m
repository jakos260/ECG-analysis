% resample.m
% function PSI=resample(PHI,tbeg,tend,nsamp)
% output samples PHI(:,tbeg:tend) resampled at nsamp equidistant points

function PSI=resampleAvo(PHI,tbeg,tend,nsamp)
dim=size(PHI);
incr=nsamp/(tend-tbeg);
x=1:tend-tbeg+1;

str=(nsamp-1)/(tend-tbeg);
x=1+(x-1)*str;
xi=1:nsamp;
PSI=interp1(x,PHI(:,tbeg:tend)',xi,'pchip');
PSI=PSI';