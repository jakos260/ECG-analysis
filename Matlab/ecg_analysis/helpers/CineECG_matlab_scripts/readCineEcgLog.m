function [DATA,RAW] = readCineEcgLog(fn)


copyfile(fn, 'c:\temp\tmp.csv');
RAW=readtable('c:\temp\tmp.csv','Delimiter',',');
delete('c:\temp\tmp.csv')
DATA= table2array( RAW(:,5:end) );    





