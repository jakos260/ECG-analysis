function colmap = mlmToqtriplotColormap( A )

colmap = A;
colmap(:,1) = (A(:,1)-1)/(length(A(:,1))-1);
