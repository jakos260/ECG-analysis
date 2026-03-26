% resample.m
% function PSI=resample(PHI,tbeg,tend,nsamp)

% output PSI is the version of PHI resampled at nsamp equidistant points

% nsamp=size(PHI,2)*sps_PSI/sps_PHI   (sps: samples per seconds)
% AvO; 2018-08-24

function PSI=resample(PHI,tbeg,tend,nsamp)

x=1:tend-tbeg+1;
str=(nsamp-1)/(tend-tbeg);
x=1+(x-1)*str;
xi=1:nsamp;
PSI=interp1(x,PHI(:,tbeg:tend)',xi,'cubic');  % 'cubic' == 'pchip'
PSI=PSI';