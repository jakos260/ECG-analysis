function values = quartile(dataIn,quart,varargin)


if size(dataIn,1)==1
    data = dataIn';
else
    data = dataIn;
end
if length(varargin) ==1
    removeVal = varargin{1};
else
    removeVal =[];
end
values = zeros(1,size(data,2));
for i=1:size(data,2)
    ds = sort(data(:,i));
    if ~isempty(removeVal)
        ds(ds==removeVal)=[];
    end
    if ~isempty(ds)
        i1= max(1,floor(quart*length(ds)));
        i2= ceil(quart*length(ds));
        di = quart*length(ds) - i1;
        values(i) = ds(i1) + (ds(i2)- ds(i1)) *di;    
    end
%     values(i) = (ds(i1)  + ds(i2))/2;    
end

