% function [cathpos_file, solflag] = cathpos(sub)
%
% This function is used to produce the txt file name of the catheter
% positions.
%
% INPUT:
% sub           = Subject number.
%
% OUTPUT:
% cathpos_file  = name txt file catheter positions.
% solflag       = [1 = true] flag to identify if readcatheterpos.m has to be run.
%
% Version 1: 01-apr-2015

function [cathpos_file, solflag] = cathpos(sub)

switch sub,
    case {7,71,72,73},
        cathpos_file    = 'cathposPig07.txt';
        solflag         = 1;
    case {8},
        cathpos_file    = 'cathposPig08.txt';
        solflag         = 1;
    case {9},
        cathpos_file    = 'cathposPig09.txt';
        solflag         = 1;
    case {10},
        cathpos_file    = 'cathposPig10.txt';
        solflag         = 1;
    case {'b01'},
        cathpos_file    = 'cathposBucket01.txt';
        solflag         = 1;
    otherwise
        cathpos_file    = '';
        warning('no catheter position txt file loaded');
        solflag         = 0;
end