function RD = RelativeDistance(A, B)
    % l = min(size(A), size(B));
    RD = norm(A - B, 'fro')/norm(A, 'fro');
end