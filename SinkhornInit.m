function [EMD,FM,TM]=SinkhornInit(A,B,WA,WB,TM)
% Data Matrices: A, B
% Weights (Disributions): WA, WB
% TM: the initialization of rotation matrix
M=size(A,2);% the number of points in set A
N=size(B,2);% the number of points in set B
% Normalize each column vector in matrices A and B
A=norCol(A);
B=norCol(B);
% Make two weight vectors be distributions
WA=WA/sum(WA);
WB=WB/sum(WB);
%% Transition Vector
meanA=sum(bsxfun(@times,WA,A),2);
meanB=sum(bsxfun(@times,WB,B),2);

A=bsxfun(@minus,A,meanA);
B=bsxfun(@minus,B,meanB);
%% Threshold
itThred=1e-4;
%% Iteration Algorithm
TM_A=norCol(TM*A);
dim=size(A,1);% the dimensionality of each point
Thre=-1;
while 1
    % Solving Eq. (9)
    % EMD is the earth move distance; FM is the flow matrix
    [EMD,FM]=Sinkhorn(TM_A,B,WA,WB);
    if abs(Thre-EMD)<=itThred
        break;
    end
    Thre=EMD;
    % Solving Eq. (10)
    resM=zeros(dim);
    temp=zeros(N,dim);
    for i=1:N
        temp(i,:)=B(1,i)*FM(:,i)'*A';
    end
    resM(1,:)=sum(temp,1);
    for i=2:dim
        resM(i,:)=sum(bsxfun(@times,(B(i,:)./B(1,:))',temp),1);
    end
    [U,~,V]=svd(resM);
    
    %resM=zeros(dim);
    %for i=1:N
    %    for j=1:M
    %        resM=resM+FM(j,i)*B(:,i)*A(:,j)';
    %    end
    %end
    %[U,~,V]=svd(resM);
    
    %MA=repmat(A,[1,N]);
    %MB=reshape(repmat(B,[M,1]),dim,[]);
    %f_ij=sqrt(FM(:)');
    %MA=bsxfun(@times,f_ij,MA);
    %MB=bsxfun(@times,f_ij,MB);
    %[U,~,V]=svd(MB*MA');
    TM=U*V';
    TM_A=TM*A;
end