function CROSS=crosscor(A,B,shift)
% crosscorrelation of A and B;
% if shift=0 :shift to zero mean is suppressed;
% dima=size(A);
% dimb=size(B);
CROSS = zeros(size(A,1),size(B,1));
for i=1:size(A,1)
    if shift~=0
        m=mean(A(i,:));
        A(i,:)=A(i,:)-m*shift;
    end
    n=norm(A(i,:));
    if n~=0
        A(i,:)=A(i,:)/n;
    end
    
    for j=1:size(B,1)
        if shift ~=0
            m=mean(B(j,:));
            B(j,:)=B(j,:)-m*shift;
        end
        n=norm(B(j,:));
        if n~=0
            B(j,:)=B(j,:)/n;
        end
        CROSS(i,j)=A(i,:)*B(j,:)';
    end
end
