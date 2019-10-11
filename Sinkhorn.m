function [EMD,FM]=Sinkhorn(TM_A,B,WA,WB)
% Squared Euclidean distance matrix
disMat=distance(TM_A,B);
disMat(disMat<0)=0;
% disMat=sqrt(disMat);
lambda=60/median(disMat(:));
K=exp(-lambda*disMat);
U=K.*disMat;
[EMD,~,u,v]=Transport(WA',WB',K,U,lambda);
FM=bsxfun(@times,v',bsxfun(@times,u,K));