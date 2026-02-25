function [ECG,names]= xmlfrompdf(fn)

X = parseXML(fn);
Y=X(2).Children(10).Children(2).Children; % children of the series
ECG=[];
names=[];
nLeads=0;
for i = 1:length(Y)    
    if strcmp(Y(i).Name,'component')
        SEQSET = Y(i).Children(2).Children;
        for j=1:length(SEQSET )
            if strcmp(SEQSET(j).Name,'component')
                Z = SEQSET(j).Children(2).Children;
                name = Z(2).Attributes(4).Value;
                if strcmp(name, 'I')
                    k= 1;
                elseif strcmp(name, 'II')
                    k= 2;
                elseif strcmp(name, 'III')
                    k= 3;
                elseif strcmp(name, 'aVR')
                    k= 4;
                elseif strcmp(name, 'aVL')
                    k= 5;
                elseif strcmp(name, 'aVF')
                    k= 6;
                elseif strcmp(name, 'V1')
                    k= 7;
                elseif strcmp(name, 'V2')
                    k= 8;
                elseif strcmp(name, 'V3')
                    k= 9;
                elseif strcmp(name, 'V4')
                    k= 10;
                elseif strcmp(name, 'V5')
                    k= 11;
                elseif strcmp(name, 'V6')
                    k= 12;
                else
                    k=-1;
%                     error('unknown lead');
                end
                if k>0
                    nLeads=nLeads+1;
                    names{k}= name;
                    scale = str2num(Z(4).Children(4).Attributes(2).Value);    
                    unit = Z(4).Children(4).Attributes(1).Value;    
                    if strcmp(unit,'uV')
                        scale =scale *1e-3;
                    end
                    ecg = str2num(Z(4).Children(6).Children.Data) * scale;
                    if isempty(ECG)
                        ECG = zeros(12,length(ecg));
                    end
                    if length(ecg) ~= size(ECG,2) 
                        warning(['number of samples per lead is not equal' num2str([size(ECG,2) length(ecg)])]);
                        ecg = [ecg(1:min(length(ecg),size(ECG,2))) (ones(1,size(ECG,2) - length(ecg))*ecg(end))];
                    end        
                    ECG(k,:) = ecg;
                end
            end
        end
    end
end

if nLeads ~= 12
    error('not all lead ECGs were stored in the xml')
end
    
%%
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end


% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end