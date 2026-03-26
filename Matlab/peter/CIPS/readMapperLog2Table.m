function TABLE = readMapperLog2Table(fn)

HEADV6948 = headerPolly6948();

if exist(fn,'file')
    status = copyfile(fn, 'c:/tmp/tmp.csv');
    TABLE  = readtable('c:/tmp/tmp.csv','Delimiter',',');
    delete('c:/tmp/tmp.csv')

    % TABLE  = readtable(fn,'Delimiter',',');
	if size(TABLE,2) == size(HEADV6948,2)
        TABLE.Properties.VariableNames = HEADV6948;
    else
       error('no matching header definition')
    end    
else
    error(['file: ' fn ' does not exist']);
end
%%
function HEADV6948 = header6948()

TABLE  = readtable('logPollyMapperheader version 6948.xlsx');
HEADV6948 = TABLE.Properties.VariableNames(2:end);