% function [ip] = inversepath(p)
%
% Depending on the input, a certain pathname is provided.
%
% INPUT & OUTPUT:
% if p = 1: ip = '/Users/arnojanssen/Documents/stw/BEM/';
% if p = 2: ip = '/Users/arnojanssen/Documents/stw/';
%
% Version 1: 01-apr-2015

function [ip] = inversepath(p)

if p == 1,
    ip = '/Users/arnojanssen/Documents/stw/PVCs/';                  % main analyses path: prepare.m & beminverse.m
elseif p == 2,
    ip = '/Users/arnojanssen/Documents/stw/';                       % path
elseif p == 3,
    ip = '/Users/arnojanssen/Documents/STW/BEM/Data/Geometries/';   % path geometries: beminverse.m
else
    error('input needed');
end