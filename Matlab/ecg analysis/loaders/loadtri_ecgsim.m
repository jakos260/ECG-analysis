function [V, F, Vidx] = loadtri_ecgsim(filepath)
%loadtri_ecgsim Read ECGSIM triangulated geometry file
%
%   [V, F, Vidx] = loadtri_ecgsim(filepath)
%
% INPUT:
%   filepath : path to geometry file
%
% OUTPUT:
%   V    : [npnt x 3] array of vertex coordinates (x,y,z)
%   F    : [ntri x 3] array of triangle indices (1-based)
%   Vidx : [npnt x 1] array of original vertex indices (optional)

    fid = fopen(filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end
    cleanup = onCleanup(@() fclose(fid));

    % ---------- Read number of points ----------
    line = nextNonEmptyLine(fid);
    npnt = sscanf(line, '%d', 1);
    if isempty(npnt) || npnt <= 0
        error('Invalid number of points in geometry file.');
    end

    V = zeros(npnt, 3);
    % Inicjalizacja Vidx tylko jeśli użytkownik go oczekuje
    if nargout > 2
        Vidx = zeros(npnt, 1);
    end

    % ---------- Read points ----------
    for i = 1:npnt
        line = nextNonEmptyLine(fid);
        vals = sscanf(line, '%f');
        if numel(vals) < 4
            error('Invalid vertex line at point %d.', i);
        end
        
        % Zapisujemy indeks wierzchołka, jeśli potrzebny
        if nargout > 2
            Vidx(i) = vals(1);
        end
        V(i, :) = vals(2:4).';
    end

    % ---------- Read number of triangles ----------
    line = nextNonEmptyLine(fid);
    ntri = sscanf(line, '%d', 1);
    if isempty(ntri) || ntri <= 0
        error('Invalid number of triangles in geometry file.');
    end

    F = zeros(ntri, 3);
    % ---------- Read triangles ----------
    for i = 1:ntri
        line = nextNonEmptyLine(fid);
        vals = sscanf(line, '%d');
        if numel(vals) < 4
            error('Invalid triangle line at triangle %d.', i);
        end
        % vals(1) is triangle index, ignore it
        F(i, :) = vals(2:4).';
    end
end

% ========================================================================
function line = nextNonEmptyLine(fid)
% Read next non-empty, non-whitespace line
    line = '';
    while isempty(line)
        if feof(fid)
            error('Unexpected end of file.');
        end
        line = strtrim(fgetl(fid));
    end
end