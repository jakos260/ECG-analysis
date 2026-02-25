% make_mla.m
% create a rectangular grid to be used by sigplot.m

% A. van Oosterom 2016_10_20


n_rows=10;

n_columns=15; 
GRID=[]; % size [rows, columns]

for i=1:n_columns,
    GRID=[GRID; [  i*ones(n_rows,1)  (1:n_rows)'       (1:n_rows)' + n_rows*(i-1) ]]
end
  
GRID=[ n_columns n_rows 0; GRID;];

 form ='%4d %4d %4d\n';
 saveasciform(['lay_' num2str(n_rows) '_' num2str(n_columns) '.mla'],GRID,form)