function Y=norCol(X)
% Normalize each column vector in matrix X
colNor=sqrt(sum(X.^2));
Y=bsxfun(@rdivide,X,colNor);