function qtriplotplay(filepath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    qtriplotplay
% 
% Desc:    Re-play a list of qtriplot commands recorded and stored by qtriplot.m 
% 
% Inputs:  filepath: path to .mat file containing variable qtriplotcmds. User is presented
%                   With a file dialog if filepath is missing or empty.
% 
% Outputs: none
% 
% Created: 05-Dec-2013 10:21:54
% 
% Author:  Peter Oosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('filepath','var') || isempty(filepath)
    [FileName,PathName] = uigetfile('*.mat','Select file containing qtriplot commands');
    filepath=fullfile(PathName,FileName);
end

load(filepath,'qtriplotcmds');
% qtriplot(qtriplotcmds);
% for n=1:40
for i=1:length(qtriplotcmds)
    if length(qtriplotcmds{i})>1 && ischar(qtriplotcmds{i}{1})
        qtriplot(sprintf(qtriplotcmds{i}{:}));
    else
        qtriplot(qtriplotcmds{i}{:});
    end
    %             pause(0.1);
end
% end