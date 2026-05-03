% file find_replicated_nodes.m
% A. van Oosterom; date: 20131011

function REPS=find_replicated_nodes(VER,small)

REPS=[];
nver=size(VER,1);

for i=1:nver;
    
    for j=i+1:nver;
        if norm(VER(i,1:3)-VER(j,1:3))<=small,
            REPS=[REPS;[i j]];
        end
    end
end
            
   
      
    
    
    