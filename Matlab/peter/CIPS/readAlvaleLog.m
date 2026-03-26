function [DATA,RAW] = readAlvaleLog(fn)


copyfile(fn, 'c:\temp\tmp.csv');
RAW = readtable('c:\temp\tmp.csv','Delimiter',',');


if size(RAW,2)==88
    DATA= table2array(RAW(:,[6:22 25:end]));   
    RAW = RAW(:,[1:22 25:end]);
else
    DATA= table2array(RAW(:,6:end));    
end

delete('c:\temp\tmp.csv')



