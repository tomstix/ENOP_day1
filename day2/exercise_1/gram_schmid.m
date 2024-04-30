function Y=gram_schmid(X) % check for linear independence
if det(X*X')==0 % input vectors are linearly dependent
    Y=[ ];
    return;
end
N=size(X,1);
Y(1,:)=X(1,:)/norm(X(1,:)); % first output vector is just a normalized first input vector
for j=2:N
    Y(j,:)=X(j,:); % initialize vector
    for k=1:j-1
        Y(j,:)=Y(j,:)-Y(k,:)*Y(j,:)'*Y(k,:); % subtract projections onto previous output vectors
    end
    Y(j,:)=Y(j,:)/norm(Y(j,:)); % normalize vector
end