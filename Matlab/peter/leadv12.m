function leadv12(varargin) %(t,PHI)

% plot standard leads

VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% cols=['b  ';'k  ';'r  ';'g  ';'c  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
cols={'b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r-.';'k-.';'g-.';'c-.';'m-.';'y-.';'r: ';'k: ';'g: ';'c: ';'m: ';'y: '};
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];
% cols = [[0.2 0.2 0.2];[0.7 0.7 0.7];[0.0,0.6,0.7]];

% cols=['k  ';'k: ';'k--';'k-.';];


warning off
paperspeed=25;
do9=0;
dowct =0;
tmax=0;
marks=0;
linew=1;
info=[];
maxphi=0;
minphi=0;
sampT=1/1000;
mycolor=[];
bw=0;
nmap=0;
Amplification=1;
PHI=[];
holdon = 0;
if length(varargin) < 1
    error('This routine needs at least two parameters');
else
    pp=1;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case '9leads'
                    do9=varargin{pp+1};pp=pp+2;
                case 'linewidth'
                    linew=varargin{pp+1};pp=pp+2;

                case 'holdon'
                    holdon=varargin{pp+1};pp=pp+2;
                case 'marks'
                    marks=varargin{pp+1};pp=pp+2;
                case 'do9'
                    do9=varargin{pp+1};pp=pp+2;
                case 'dowct'
                    dowct=varargin{pp+1};pp=pp+2;
                case 'info'
                    info=varargin{pp+1};pp=pp+2;
                case 'max'
                    maxphi=varargin{pp+1};pp=pp+2;
                case 'bw'
                    bw=varargin{pp+1};pp=pp+2;

                case 'sampt'
                    sampT=varargin{pp+1};pp=pp+2;
                case 'paperspeed'
                    paperspeed=varargin{pp+1};pp=pp+2;
                case 'amplification'
                    Amplification=varargin{pp+1};pp=pp+2;
                case 'color'
                    mycolor = varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        else
            if iscell(varargin{pp})

                for i=1:length(varargin{pp})
                    nmap = nmap+1;
                    PHI{nmap} = varargin{pp}{i};
                    tmax = max(tmax,size(PHI{pp},2)-1);
                end
            else
                nmap = nmap+1;

                PHI{nmap} = varargin{pp};
                tmax = max(tmax,size(PHI{pp},2)-1);

            end

            pp=pp+1;
        end
    end
end
sampT=1000*sampT;
tmax=tmax*sampT+1*sampT;

if length(maxphi)>1
    minphi=maxphi(1);
    maxphi=maxphi(2);
elseif maxphi~=0
    minphi=-maxphi;
end


%% wct reference and determine maxphi
wct = 1:3;
vs =  [wct 4:9];

for i=1:nmap
    A=PHI{i};
    if dowct
        A = A - ones(size(A,1),1)*mean(A(wct,:));
    end
    PHI{i} = A;
    if maxphi==0
        maxphi=max(maxphi,max(PHI{i}(:)));
        minphi=min(minphi,min(PHI{i}(:)));
    end
end

maxphi = max(abs([maxphi minphi]));
maxphi = ceil(maxphi * 2)/2;
minphi = -maxphi;

if paperspeed<=50
    tmax=ceil(tmax/paperspeed)*paperspeed;
end

if ~holdon, clf; end

fs=11;

fweight='demi';


if do9
    ny=3;nx=3;
    x=0;
else
    x=0;
    ny=3;nx=4;
    textI=['I  ';'II ';'III'];
    for k=1:3
        y=3-icyc(k,3);
        axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
        for i=1:nmap%:-1:1
            hold on
            if size(PHI{i},1) == 12
                lead1 = PHI{i}(1,:);
                lead2 = PHI{i}(2,:);
                lead3 = PHI{i}(3,:);
            else
                phi1=PHI{i}(wct(1),:);
                phi2=PHI{i}(wct(2),:);
                phi3=PHI{i}(wct(3),:);
                lead1 = phi2-phi1;
                lead2 = phi3-phi1;
                lead3 = phi3-phi2;
            end
            if ~isempty(mycolor)
                if size(mycolor,1)==1
                    color = mycolor(:);
                else
                    color = mycolor(i,:);
                end
            else
                color = cols{i};
            end

            if ~bw
                if k==1
                    plot(sampT*(0:length(lead1)-1), lead1,color,'Linewidth',linew);
                elseif k==2
                    plot(sampT*(0:length(lead2)-1), lead2,color,'Linewidth',linew);
                else
                    plot(sampT*(0:length(lead3)-1), lead3,color,'Linewidth',linew);
                end
            else
                if k==1
                    plot(sampT*(0:length(lead1)-1), lead1,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                elseif k==2
                    plot(sampT*(0:length(lead2)-1), lead2,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                else
                    plot(sampT*(0:length(lead3)-1), lead3,styles(icyc(i,length(styles)),:),'Linewidth',linew);
                end
            end
        end
        axis([0,tmax,minphi,maxphi]); axis off
        if ~holdon
            text(tmax*0.65,0.9*maxphi,textI(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
        end
        ecgraster('SampleRate',1000/sampT,'Paperspeed',paperspeed,'Amplification',Amplification,'drawgrid',~holdon)

    end
    x=1;
end

for k=1:9
    y=3-icyc(k,3);
    axes('Position',[x/nx+0.005,y/ny+.03,1/nx-0.015,.9/ny]);
    if mod(k,ny)==0, x=x+1; end
    for i=1:nmap%:-1:1
        doCorrect =1;
        if size(PHI{i},1) == 12
            doCorrect  = 0;
            PHI9 = PHI{i}(4:end,:);
        else
            PHI9 = PHI{i};
        end
        V=vs(k);
        phi=PHI9(V,:);
        if k >= 1 && k <=3 && doCorrect
            phi= phi* 1.5;
        end
        if ~isempty(mycolor)
            if size(mycolor,1)==1
                color = mycolor(:);
            else
                color = mycolor(i,:);
            end
        else
            color = cols{i};
        end
        if ~bw
            plot(sampT*(0:length(phi)-1), phi,color,'Linewidth',linew);
        else
            plot(sampT*(0:length(phi)-1), phi,styles(icyc(i,length(styles)),:),'Linewidth',linew);
        end
        hold on;
    end
    axis([0,tmax,minphi,maxphi]); axis off
    if ~holdon
        text(tmax*0.5,0.9*maxphi,VLab(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
    end
    ecgraster('Amplification',Amplification,'SampleRate',1000/sampT,'Paperspeed',paperspeed,'drawgrid',~holdon)
end
annotation('textbox',[0.65 0.01 0.6 0.025],'string',[num2str(paperspeed) ' mm/s  ' num2str(10/Amplification) ' mm/mV'],'edgecolor','none','FitBoxToText','on');
end
%%
function ecgraster(varargin)
Amplification=1; % 1mV/mm
Paperspeed=25; % 25 mm/s
Calibrate='off';
Scale=1;
UseAxis=gca;
SampleRate=1000; %
drawgrid =1;
ProcessVararginHere(varargin);
Paperspeed=Paperspeed/SampleRate;

col =[1,.5,.5];
col = [132, 232, 177]/255;


axes(UseAxis);
c0=get(UseAxis,'children');
hold on
if strcmpi(Calibrate,'on')
    % oldunits=get(gca,'Units');
    set(gca,'Units','centimeters');
    yl=ylim;xl=xlim;
    xw=diff(xlim)*Paperspeed/10*Scale;
    yh=diff(ylim)/Amplification*Scale;
    pos=get(gca,'position');
    set(gca,'xlim',xl,'ylim',yl,'Position',[pos(1),pos(2),xw,yh]);
    set(gca,'Units',oldunits);
end
daspect([5/Paperspeed,0.5*Amplification,1]);

yl=ylim;
xl=xlim;
if drawgrid
    if diff(yl)/(.1*Amplification)<1000
        for k=yl(1):.1*Amplification:yl(2);
            plot(xl,k*[1,1],'color',col,'linewidth',.5)
        end
        for k=yl(1):.5*Amplification:yl(2);
            plot(xl,k*[1,1],'color',col,'linewidth',1.5)
        end
    else
        disp(['Too many lines in ',mfilename,', possibly wrong amplitude']);
    end
    if diff(xl)*Paperspeed<1000
        for k=xl(1):1/Paperspeed:xl(2);
            plot(k*[1,1],yl,'color',col,'linewidth',.5)
        end
        %    for k=xl(1):5/Paperspeed:xl(2);
        %       plot(k*[1,1],yl,'color',[1,.5,.5],'linewidth',1)
        %    end
        for k=xl(1):5/Paperspeed:xl(2);
            plot(k*[1,1],yl,'color',col,'linewidth',1.5)
        end

    else
        disp(['Too many lines in ',mfilename,', possibly wrong sample freq']);
    end
end
set(UseAxis,'children',[c0;setdiff(get(UseAxis,'children'),c0)]);
end
%%

function ProcessVararginHere(args,varargin)
if isempty(args);return;end
Hide = [];
ProcessVararginHere(varargin);
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
end
%%
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
end

%%
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

end
%%
% function i=icyc(k,n)
% k and n are integers;  n a scalar
% 20110204
% rem type of function; example:
% for n=3 and
% k=[ -3 -2 -1 0 1 2 3 4 5 6 7 etc ]; it produces
% i=[  3  1  2 3 1 2 3 1 2 3 1 etc;;
function k=icyc(k,n)
k=rem(k,n);
k(k<=0)= k(k<=0)+n;
end


%%

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
        % input is a java xml object
        xDoc = file;
    else
        %check for existance
        if (exist(file,'file') == 0)
            %Perhaps the xml extension was omitted from the file name. Add the
            %extension and try again.
%             if (isempty(strfind(file,'.xml')))
%                 file = [file '.xml'];
%             end
            
            if (exist(file,'file') == 0)
                error(['The file ' file ' could not be found']);
            end
        end
        %read the xml file
        xDoc = xmlread(file);
    end
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    
end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = struct; textflag = 'Text';
    if hasChildNodes(theNode)
        childNodes = getChildNodes(theNode);
        numChildNodes = getLength(childNodes);

        for count = 1:numChildNodes
            theChild = item(childNodes,count-1);
            [text,name,attr,childs,textflag] = getNodeData(theChild);
            
            if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
                    if (~iscell(children.(name)))
                        %put existsing element into cell format
                        children.(name) = {children.(name)};
                    end
                    index = length(children.(name))+1;
                    %add new element
                    children.(name){index} = childs;
                    if(~isempty(fieldnames(text)))
                        children.(name){index} = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name){index}.('Attributes') = attr; 
                    end
                else
                    %add previously unknown (new) element to the structure
                    children.(name) = childs;
                    if(~isempty(text) && ~isempty(fieldnames(text)))
                        children.(name) = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name).('Attributes') = attr; 
                    end
                end
            else
                ptextflag = 'Text';
                if (strcmp(name, '#cdata_dash_section'))
                    ptextflag = 'CDATA';
                elseif (strcmp(name, '#comment'))
                    ptextflag = 'Comment';
                end
                
                %this is the text in an element (i.e., the parentNode) 
                if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                    if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                        ptext.(ptextflag) = text.(textflag);
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                    end
                end
            end
            
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = toCharArray(getNodeName(theNode))';
    name = strrep(name, '-', '_dash_');
    name = strrep(name, ':', '_colon_');
    name = strrep(name, '.', '_dot_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text,textflag] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)) && isempty(fieldnames(text)) && ~isempty(getTextContent(theNode)))
        %get the data of any childless nodes
        % faster than if any(strcmp(methods(theNode), 'getData'))
        % no need to try-catch (?)
        % faster than text = char(getData(theNode));
        text.(textflag) = toCharArray(getTextContent(theNode))';
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if hasAttributes(theNode)
       theAttributes = getAttributes(theNode);
       numAttributes = getLength(theAttributes);

       for count = 1:numAttributes
            %attrib = item(theAttributes,count-1);
            %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
            %attributes.(attr_name) = char(getValue(attrib));

            %Suggestion of Adrian Wanner
            str = toCharArray(toString(item(theAttributes,count-1)))';
            k = strfind(str,'='); 
            attr_name = str(1:(k(1)-1));
            attr_name = strrep(attr_name, '-', '_dash_');
            attr_name = strrep(attr_name, ':', '_colon_');
            attr_name = strrep(attr_name, '.', '_dot_');
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end