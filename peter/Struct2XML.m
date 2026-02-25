 function struct2xml(filename,structname) 
global i_struct2xml 
fid = fopen(filename,'w'); 
disp(filename); 
fprintf(fid,['<?xml version=',34,'1.0',34,' ?>\n']); 
evalin('caller','global i_struct2xml'); 
i_struct2xml.level = -1; 
evalin('caller',sprintf('i_struct2xml.struct=%s;',structname)); 
xmlprint(fid,structname,i_struct2xml.struct); 
evalin('caller','clear global i_struct2xml'); 
fclose(fid); 
%----------------------- 
function xmlprint(fid,structname,struct) 
global i_struct2xml 
i_struct2xml.level = i_struct2xml.level+1; 
if i_struct2xml.level==0 
   spaces = ''; 
else 
   spaces = char(ones(i_struct2xml.level,1) * 9); 
end 

if isstruct(struct) 
	if length(struct)>1
		i_struct2xml.level = i_struct2xml.level-1;
		for k=1:length(struct)
			xmlprint(fid,structname,struct(k));
		end
		return
	end
   if isfield(struct,'Attributes') 
      fprintf(fid,'%s<%s',spaces,structname); 
      flds = fieldnames(getfield(struct,'Attributes')); 
      for k = 1:length(flds) 
         fprintf(fid,[' %s=',34,'%s',34,''],flds{k},getfield(getfield(struct,'Attributes'),flds{k})); 
      end 
      if isfield(struct,'Value') 
         fprintf(fid,'>%s</%s>\n',struct.Value,structname); 
         i_struct2xml.level = i_struct2xml.level-1; 
         return 
      else 
         fprintf(fid,'>\n'); 
      end 
 
   else 
      fprintf(fid,'%s<%s>\n',spaces,structname); 
   end 
 
   flds = fieldnames(struct); 
   for k = 1:length(flds) 
      if ~strcmp(flds(k),'Attributes') 
         if ~(isfield(struct,'Attributes') & strcmp(flds(k),'Value')) 
            xmlprint(fid,flds{k},getfield(struct,flds{k})); 
         end 
      end 
   end 
   fprintf(fid,'%s</%s>\n',spaces,structname); 
elseif size(struct,1)==0 
   fprintf(fid,'%s<%s></%s>\n',spaces,structname,structname); 
elseif iscell(struct) 
   l = size(struct,2); 
   i_struct2xml.level = i_struct2xml.level-1; 
   for k = 1:l 
      cl = struct{k}; 
      xmlprint(fid,structname,cl); 
   end 
   i_struct2xml.level = i_struct2xml.level+1; 
 
elseif ischar(struct) 
   l = size(struct,1); 
   for k = 1:l 
      item = struct(k,:); 
       item = strrep(item,'&','&amp;'); 
item = strrep(item,'<','&lt;'); 
item = strrep(item,'>','&gt;'); 
item = strrep(item,char(34),'&quot;'); 
item = strrep(item,char(39),'&apos;'); 
 
 
      fprintf(fid,'%s<%s>%s</%s>\n',spaces,structname,item,structname); 
   end 
else %isnum 
   l = prod(size(struct)); 
   t = sprintf('%.10g ',struct(:)'); 
   fprintf(fid,'%s<%s type=%cnum%c>%s</%s>\n',spaces,structname,34,34,t,structname); 
end 
 
i_struct2xml.level = i_struct2xml.level-1; 
 

