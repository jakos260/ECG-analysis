for i=1:size(blindlog,2)
    [pathstr, name, ext] = fileparts(blindlog(i).exportfilename);
    display(sprintf('%d:\t%s.\t(%s)',i,name,blindlog(i).comment));
end;