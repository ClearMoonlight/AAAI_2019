clear
clc

d=500;
n1=20000;
n2=30000;
tarD=50;% the dimensionality of the manifold
%% Low-Dimensional Manifold
disp('two different manifolds');
W=randn(d);
A=zeros(d,n1);
B=zeros(d,n2);
pA=randn(tarD,n1);
pB=randn(tarD,n2);
R1=randi(d,1,tarD);
R2=randi(d,1,tarD);
A(R1,:)=pA;
B(R2,:)=pB;
A=W*A;
B=W*B;

WA=rand(1,n1);
WB=rand(1,n2);

IT=randn(d);
%% Original EMD
disp('Original EMD');
tic
[EMD1,~,~]=SinkhornInit(A,B,WA,WB,IT);
EMD1
T1=toc;
%% EMD based on Core-Set
disp('EMD based on Core-Set');
csArr=[50,40,30,20,10];% the array of core-set size
csSz=numel(csArr);% the number of core-sets
runNo=20;
emdTab=zeros(csSz,runNo);% record emd
timTab=zeros(csSz,runNo);% record time

copA=norCol(A);
copB=norCol(B);
copWA=WA/sum(WA);
copWB=WB/sum(WB);
meanA=sum(bsxfun(@times,copWA,copA),2);
meanB=sum(bsxfun(@times,copWB,copB),2);
copA=bsxfun(@minus,copA,meanA);
copB=bsxfun(@minus,copB,meanB);
for outI=1:csSz
    k1=ceil(n1/csArr(outI));% the size of core-set of set A
    k2=ceil(n2/csArr(outI));% the size of core-set of set B
    for inI=1:runNo
        outI,inI
        tic
        [CSA,CSWA]=KCenter(A,WA,k1);
        [CSB,CSWB]=KCenter(B,WB,k2);
        [~,~,TM2]=SinkhornInit(CSA,CSB,CSWA,CSWB,IT);
        [EMD2,~]=Sinkhorn(TM2*copA,copB,copWA,copWB);
        timTab(outI,inI)=toc;
        emdTab(outI,inI)=EMD2;
        EMD2
    end
end
save('Result/RD/LDM/1/EMD1.mat','EMD1');
save('Result/RD/LDM/1/T1.mat','T1');
save('Result/RD/LDM/1/emdTab.mat','emdTab');
save('Result/RD/LDM/1/timTab.mat','timTab');