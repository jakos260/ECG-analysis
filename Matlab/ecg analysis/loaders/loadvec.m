function A = loadvec(name)
    f=fopen(name);
    if (f==-1)
      fprintf('\nCannot open %s\n', name);
      A=0;
      return;
    end

    % Read matrix size
    dims = fscanf(f, '%d %d', 2);
    rows = dims(1);
    cols = dims(2);
    
    % Read the data
    data = fscanf(f, '%f');
    
    fclose(f);
    
    % Reshape into matrix
    A = reshape(data, rows, cols);

end