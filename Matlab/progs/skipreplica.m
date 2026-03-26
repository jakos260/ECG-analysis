% skipreplica.m
% SINGLES=skipreplica(LIST,column,small,mode)
% FIRST: rows of LIST are sorted based on the values in the column: column
% (default 1)
% NEXT: skipreplica skips any duplicate rows
% if mode==0, preliminary sorting is skipped; takes longer (?)
% default skipping criterium: max(abs(entries))<small

% A. van Oosterom; 2006/01/19
% 2007_04_05: see alternative script

function SINGLES=skipreplica(LIST,column,small,mode);
  if nargin==1, column=1; end
  if nargin<=2,small=eps; end 
  if nargin<=3,mode=1; end
  SINGLES=[];
  if isempty(LIST)==1, return, end
  [ni nj]=size(LIST); transpos=0;
  if ni==1, LIST=LIST';transpos=1; end
  
  
%   if mode==1, 
%      LIST=sortrows(LIST,column);
%      SINGLES(1,:)=LIST(1,:);
%      len=size(LIST,1);
%      if len==1,
%          return,
%      end
%      for i=2:len,
%         if abs(LIST(i,column)-LIST(i-1,column))>small, SINGLES=[SINGLES; LIST(i,:)];,end
%      end
%   else, 
%      while isempty(LIST)==0,
%            single=LIST(1,:);
%            SINGLES=[SINGLES;single];
%            list=find(abs(LIST(:,column)-single(column))<=small);
%            LIST(list,:)=[];
%       end
%   end

LIST=sortrows(LIST,column);

% SINGLES=LIST(1,:);
% LIST(1,:)=[];
% while isempty(LIST)==0,
%     SINGLES=[SINGLES;LIST(1,:)];
%     LIST(1,:)=[];
%     k=find(sum(LIST == ones(size(LIST,1),1)*LIST(1,:),2)==2);
%     if isempty(k)==0,
%         LIST(k,:)=[];
%     end
% end
%  
         
 if transpos==1, SINGLES=SINGLES';end
