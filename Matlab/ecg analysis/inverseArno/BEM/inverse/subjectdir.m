% function [subject, bsmdir, layfile] = subjectdir(sub, bp, issinus)
%
% This function is used to produce the directory and file names for a
% specific subject and study.
%
% INPUT
% sub:          Subject number.
% bp:           Main directory path for specific study. Defined in inversepath.m
% issinus:      ...
%
% OUTPUT
% subject:    	Name of subject [folder].
% bsmdir:   	Directory where results from prepare.m are stored [.ecgspecs & .selecg].
% layfile:      Layout file name.
%
% Version 1: 01-apr-2015

function [subject, bsmdir, layfile] = subjectdir(sub, bp, issinus)

if ~exist('issinus', 'var'), issinus = 0; end
if ~exist('bp', 'var'), warning('no bp given, bp = '';'); bp = ''; end

studyname       = 'NICE00';
dirbeats        = '/Biosemi/export/beats/';
subject         = [studyname num2str(sub)];
bsmdir          = [bp subject dirbeats];
layfile         = 'ams65.mla';