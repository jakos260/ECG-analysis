%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [orgnames,surface,stimelectrodes,actmaps] = readcatheterpos(filename)
%
% Name:    readcatheterpos
%
% Desc:    Read information on catheter positions from txt file (STWCommnoData)
%
% Inputs:  filename: Filename or path to tab deilimited input file
%
% Outputs: orgnames:        filename parts for each data block
%          surface:         surface of stimulation (EPI, LVEndo, RVEndo)
%          stimelectrodes:  position of stimulation electrode. [Exact, nearest point on surface]
%          actmaps:         catheter positions ([Exact, surface] and activation times rel. to stim
%
% Created: 07-Oct-2013 17:25:42
%
% Author:  peteroosterhoff
%
% Version 1: 01-apr-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orgnames,surface,stimelectrodes,actmaps] = readcatheterpos(filename)

str         = fileread(filename);
str         = regexprep(str,'''','');                           	% remove single quotes
str         = regexprep(str,'"','');                             	% remove double quotes
linestr     = regexp(str,'[\r\n]','split');
dellines    = ones(1,length(linestr));                              % Remove comment(end) lines

for i = 1:length(linestr),
    if linestr{i}(1) == '#' || linestr{i}(1) == '"', dellines(i) = 0; end
end

linestr = linestr(find(dellines));

ind     = 1;
lind    = 1;

while lind <= length(linestr)-2,                                    % minimum map
    orgnames{ind}   = skipempty(regexp(linestr{lind},'\t','split'));
    lind            = lind+1;
    
    surface{ind}    = regexp(linestr{lind},'(\w*){1}','match');
    lind            = lind+1;
    C               = skipempty(textscan(linestr{lind},'%f',6,'Delimiter','\t'));
    
    if length(C{1}) == 6,
        stimelectrodes(ind).exactpos    = C{1}(1:3)';
        stimelectrodes(ind).surfpos     = C{1}(4:6)';
        lind                            = lind+1;
    else
        error('parse error');
    end
    
    actmaps(ind).nmap   = 0;                                        % number op measurements in map
    stopflag            = false;
    
    while ~stopflag && lind <= length(linestr),
        C = skipempty(textscan(linestr{lind},'%s %s %f %f %f %f %f %f %d %d %d %f',1,'Delimiter','\t','MultipleDelimsAsOne',1));
        
        if length(C) ~= 12 && length(C) ~= 11,
            if isempty(C), stopflag = true; end                   	% empty line, start of next block
            
            lind = lind+1;                                          % ignore incomplete lines
            
        else                                                        % not break or error
            actmaps(ind).nmap                           = actmaps(ind).nmap+1;
            actmaps(ind).surf{actmaps(ind).nmap}        = C{1}{1};
            actmaps(ind).label{actmaps(ind).nmap}       = C{2}{1};
            actmaps(ind).exactpos(actmaps(ind).nmap,:)  = cell2mat(C(3:5));
            actmaps(ind).surfpos(actmaps(ind).nmap,:)   = cell2mat(C(6:8));
            
            % manshift not exported
            if length(C) >= 12,
                actmaps(ind).dep(actmaps(ind).nmap) = cell2mat(C(12));
            else
                actmaps(ind).dep(actmaps(ind).nmap) = NaN;          % give electrode position even is act time missing
            end
            
            lind = lind+1;
        end
        
    end
    
    ind = ind+1;
end

% Sub-function:
function cellstrout = skipempty(cellstrin)
cellstrout = {};

for i = 1:length(cellstrin),
    if ~isempty(cellstrin{i}), cellstrout{end+1} = cellstrin{i}; end
end