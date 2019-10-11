clear
clc

dim=10;
M=1000;
N=5000;
FM=rand(M,N);
A=rand(dim,M);
B=rand(dim,N);
%% first try
resM1=zeros(dim);
for i=1:N
    for j=1:M
        resM1=resM1+FM(j,i)*B(:,i)*A(:,j)';
    end
end
resM1
%% second try
resM2=zeros(dim);
temp=zeros(N,dim);
for i=1:N
    temp(i,:)=B(1,i)*FM(:,i)'*A';
end
resM2(1,:)=sum(temp,1);
for i=2:dim
    resM2(i,:)=sum(bsxfun(@times,(B(i,:)./B(1,:))',temp),1);
end
resM2
%% third try
B*FM'*A'