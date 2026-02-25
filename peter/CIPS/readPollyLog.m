function [DATA,RAW] = readPollyLog(fn)


copyfile(fn, 'tmp.csv');
RAW  = readtable('tmp.csv','Delimiter',',');
DATA = table2array(RAW(:,10:end));    

delete('tmp.csv')



