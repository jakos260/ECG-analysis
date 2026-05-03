% getpeaks.m
% function PEAKS=getpeaks(signal,threshold,window)
% window number of samples in narrow interval around maxima
% Find values and locations of local maxima(peaks) in a (scalar) signal.
% dedicated to finding RR-intervals in a record containing several beats.
% in nargin<3, all signal values above a threshold are considered
% if nargin==3, just the first value in the widows (samples) around the maximum are
% documented, if window is entered negative, it is the last vaue peak within the window
% is documented
% last modified: 2013_01_08; A. van Oosterom

function PEAKS=getpeaks(signal,threshold,window)

nt=size(signal,2);

signal(signal<threshold)=0;
PEAKS=[];
% treat possible extremes at the bounds of the interval
augmented=[ 0 signal 0 ];

SEARCH=[(1:nt+1)' augmented(1:nt+1)' augmented(2:nt+2)'];

starters=SEARCH(SEARCH(:,2)==0 & SEARCH(:,3)~=0);
endings=SEARCH(SEARCH(:,2)~=0 & SEARCH(:,3)==0)-1;
npeaks=size(starters,1);

if npeaks>0,
    for i=1:npeaks,
        [PEAKS(i,1) PEAKS(i,2)]=max(signal(starters(i):endings(i)));
    end
    PEAKS(:,2)=PEAKS(:,2)+ starters-1;
end

if nargin<3, return, end

SEARCH=[];

if sign(window)>0,
    while ~isempty(PEAKS),
        SEARCH=[SEARCH;PEAKS(1,:)];
        PEAKS(PEAKS(:,2)< PEAKS(1,2)+window,:)=[];
    end
    PEAKS=SEARCH;
else,
    while ~isempty(PEAKS),
        SEARCH=[SEARCH;PEAKS(end,:)];
        PEAKS(PEAKS(:,2)>PEAKS(1,2)-window,:)=[];
    end
    PEAKS=reverse(SEARCH);
end    
    
    
    
    
    
    
    
