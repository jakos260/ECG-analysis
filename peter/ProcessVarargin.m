 function ProcessVarargin(args,varargin) 
if isempty(args);return;end 
Hide = []; 
ProcessVarargin(varargin); 
if isempty(Hide) 
   Hide = 'varargin'; 
elseif iscell(Hide) 
   Hide  =  [Hide,{'varargin'}]; 
else 
   Hide  =  [{Hide},{'varargin'}]; 
end 
call_who = setdiff(evalin('caller','who'),Hide); 
stack = dbstack; 
k = 1; 
while k <= length(args) 
   if ~isempty(strmatch(strtok(args{k},'.'),call_who,'exact')) 
      if k==length(args) 
         error(['unpaired name/value ',args{k},' for ',stack(2).name]); 
      end 
      if any(args{k}=='.') 
         assignin('caller','cxqweryu',args{k+1}); 
         evalin('caller',[args{k},'=cxqweryu;clear cxqweryu;']); 
      else 
         assignin('caller',args{k},args{k+1}); 
      end 
      k = k+2; 
   else 
      switch (args{k}(1)) 
      case '$' 
         [p,f,e] = fileparts(args{k}(2:end)); 
         disp(f); 
         switch (e) 
         case '.mat' 
            evalin('caller',['load(''',p,''');']); 
            k  =  k+1; 
         case '.m' 
            oldcd  =  evalin('caller','cd') ; 
            evalin('caller',sprintf('cd %s;%s;cd %s',p,f,oldcd)); 
            k  =  k+1; 
         case '.xml' 
            args2  =  ParseXML(xml2struct(args{k}(2:end))); 
            args  =  [args(1:k-1),args2,args(k+1:end)]; 
         otherwise 
            disp(['unsupported file format: ',e]); 
            k  =  k+1; 
         end 
      case '!' 
         evalin('caller',args{k}(2:end)); 
         k = k+1; 
      otherwise 
         disp(['valid variables for ',stack(2).name,':']) 
         disp(call_who) 
         error(['Illegal argument :',args{k},' for ',stack(2).name]); 
      end 
   end 
end 
 
function args  =  ParseXML(xmlstruct) 
args  =  []; 
flds  =  fieldnames(xmlstruct); 
for k  =  1:length(flds) 
   fld  =  getfield(xmlstruct,flds{k}); 
   if strcmp(lower(flds{k}),'options') 
      args  =  [args,ParseXMLOptions(fld)]; 
   else 
      if isstruct(fld) 
         args  =  [args,ParseXML(fld)]; 
      end 
   end 
end 
 
function args  =  ParseXMLOptions(xmlstruct) 
args  =  []; 
flds  =  fieldnames(xmlstruct); 
for k  =  1:length(flds) 
   fld  =  getfield(xmlstruct,flds{k}); 
   if iscell(fld) | isstruct(fld) 
      if iscell(fld) 
         args2 = []; 
         for k1 = 1:length(fld) 
            args2(k1 * 2+[-1,0])  =  ParseXMLOptions(fld{k1}); 
         end 
      else 
         args2  =  ParseXMLOptions(fld); 
      end 
      for k1  =  1:2:length(args2) 
         args  =  [args,{[flds{k},'.',args2{k1}],args2{k1+1}}]; 
      end 
   else 
      args  =  [args,{flds{k},fld}] ; 
   end 
end 
 
