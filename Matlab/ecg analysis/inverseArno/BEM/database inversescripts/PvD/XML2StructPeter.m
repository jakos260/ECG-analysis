function OUTPUT = XML2StructPeter(filename) 
global i_xml2struct 
if ~strcmp(lower(filename(end-3:end)),'.xml') 
   filename = [filename,'.xml']; 
end 

fid = fopen(filename,'rt'); 
txt = char(fread(fid)'); 
fclose(fid); 

xmlheader = ['<?xml version=',34,'1.0',34,' ']; %char(34)=double quote 
xml0 = max([1,findstr(txt,xmlheader)+length(xmlheader)]); 
if xml0==1 
   xmlheader = '<?xml version = ''1.0'' ?>'; 
   xml0 = max([1,findstr(txt,xmlheader)+length(xmlheader)]); 
   if xml0==1 
      disp('not a valid xml file'); 
      return 
   else 
      %assume this is a mtangled xml file with only single quotes 
      txt = strrep(txt,'''',char(34)); 
   end 
end 
i_xml2struct.struct = readXML(txt(xml0:end)); 
if nargout==0 
   name = fieldnames(i_xml2struct.struct);name = name{1}; 
   evalin('caller',sprintf('global i_xml2struct;%s=i_xml2struct.struct.%s;clear i_xml2struct;',name,name)); 
   OUTPUT = name; 
else 
   OUTPUT = i_xml2struct.struct; 
end 



%----------------------- 
function [OUTPUT,last] = readXML(txt) 
global i_xml2struct 
OUTPUT = []; 
first = 1; 
txtlen = length(txt); 
last = 0; 
while (1) 
   i = findstr(txt,'<'); 
   if isempty(i) 
      OUTPUT = txt; 
      last = length(txt); 
      return 
   end 
   i = i(1); 
   tag = strtok(txt(i+1:end),'>'); 
   first = i+length(tag)+2; 
   
   if any(isspace(tag)) 
      [tag,attributes] = parse_attributes(tag); 
   else 
      attributes = []; 
   end 
   
   iopen = sort([findstr(txt,['<',tag,' ']),findstr(txt,['<',tag,'>'])]); 
   iclose = findstr(txt,['</',tag,'>']); 
   if length(iopen)~=length(iclose) 
      if ~isempty(findstr(txt,['/>'])) %hmm, may not always work 
         iclose = iopen; 
         tag = tag(1:end-1); 
      else 
         error(['open and close tags do not match for:',tag]); 
         return 
      end 
   end 
   s = [ones(size(iopen)),-ones(size(iclose))]; 
   [ioc,idx] = sort([iopen,iclose]); 
   s = cumsum(s(idx)); 
   last = ioc(min(find(s==0)))-1; 
   [item,last] = readXML(txt(first:last)); 
   if ischar(item) 
      item = strrep(item,'&lt;','<'); 
      item = strrep(item,'&gt;','>'); 
      item = strrep(item,'&quot;',char(34)); 
      item = strrep(item,'&apos;',char(39)); 
      item = strrep(item,'&amp;','&'); 
      
      
      idx = find(item~=32 & item~=10 & item ~=13 & item~=9); 
      item = item(min(idx):max(idx)); 
      
      
   end 
   
   if ~isempty(attributes) 
      if length(fieldnames(attributes))==1 & isfield(attributes,'type') 
         if strcmp(attributes.type,'num') 
            item(find(item==10)) = 32; 
            item(find(item==13)) = 32; 
            item = str2num(item); 
         end 
      else 
         if ischar(item);it=item;item=[];item.Value  =  it;end 
         item  =  setfield(item,'Attributes',attributes); 
      end 
   end 
   
   last = last+first-1; 
   if isfield(OUTPUT,tag) 
      tmp = getfield(OUTPUT,tag); 
      if ischar(tmp) 
         l = size(tmp,1); 
      else 
         l = length(tmp); 
      end 
      
      if (l==1) 
         OUTPUT = setfield(OUTPUT,tag,{tmp,item}); 
      else 
         eval(sprintf('OUTPUT.%s{%d}=item;',tag,l+1)); 
      end 
      
   else 
      if findstr(tag,'/ke')
         keyboard
      end
      if  ~isempty(strfind(tag,'.'))
          tag(strfind(tag,'.'))='_'; %remove invalid filename tokens
      end
      OUTPUT = setfield(OUTPUT,tag,item); 
   end 
   txt = txt(last+length(tag)+4:end); 
   if isempty(txt);break;end; 
   if all(isspace(txt));break;end; 
end 
last = txtlen; 


%---------------------- 
function [tag,attrib] = parse_attributes(tagplus)
[tag,atttxt] = strtok(tagplus); 
attrib = []; 
while ~isempty(atttxt) 
   [attname,atttxt] = strtok(atttxt,'='); 
   atttxt = atttxt(2:end); % skip = 
   attname = strtok(attname); 
   atttxt = fliplr(deblank(fliplr(atttxt))); 
   atttxt = atttxt(min([1,find(isspace(atttxt))]):end); 
   [attval,atttxt] = strtok(atttxt,atttxt(1)); 
   atttxt = atttxt(2:end); %skip ' or ' 
   attname(strfind(attname,':'))='_'; %remove invalid filename tokens
   attname(strfind(attname,'/'))='_'; %remove invalid filename tokens
   attname(strfind(attname,'\'))='_'; %remove invalid filename tokens
   attrib = setfield(attrib,attname,attval); 
end 


