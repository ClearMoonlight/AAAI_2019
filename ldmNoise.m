clear
clc

load('data/A.mat');
load('data/B.mat');
load('data/WA.mat');
load('data/WB.mat');

dim=size(A,1);% the dimensionality of point
n1=size(A,2);% the number of points in set A
n2=size(B,2);% the number of points in set B
diaA=sqrt(max(max(distance(A))));% the diameter of set A
diaB=sqrt(max(max(distance(B))));% the diameter of set B

csArr=[50,40,30,20,10];% the array of core-set size
csSz=numel(csArr);% the number of core-sets
ratArr=[0.005,0.010,0.015,0.020,0.025];%[0.03,0.06,0.09,0.12,0.15];
stA=diaA*ratArr;
stB=diaB*ratArr;
noiNum=numel(ratArr);
runNo=20;

IT=randn(dim);
for i=1:noiNum
    i
    Ai=A+stA(i)*randn(dim,n1);
    Bi=B+stB(i)*randn(dim,n2);
    %% Original data
    tic
    [EMD1,~,~]=SinkhornInit(Ai,Bi,WA,WB,IT);
    T1=toc;
    %% Core-sets
    copA=norCol(Ai);
    copB=norCol(Bi);
    copWA=WA/sum(WA);
    copWB=WB/sum(WB);
    meanA=sum(bsxfun(@times,copWA,copA),2);
    meanB=sum(bsxfun(@times,copWB,copB),2);
    copA=bsxfun(@minus,copA,meanA);
    copB=bsxfun(@minus,copB,meanB);
    emdTab=zeros(csSz,runNo);
    timTab=zeros(csSz,runNo);
    for outI=1:csSz
        k1=ceil(n1/csArr(outI));% the size of core-set of set A
        k2=ceil(n2/csArr(outI));% the size of core-set of set B
        for inI=1:runNo
            outI,inI
            tic
            [CSA,CSWA]=KCenter(Ai,WA,k1);
            [CSB,CSWB]=KCenter(Bi,WB,k2);
            [~,~,TM2]=SinkhornInit(CSA,CSB,CSWA,CSWB,IT);
            [EMD2,~]=Sinkhorn(TM2*copA,copB,copWA,copWB);
            timTab(outI,inI)=toc;
            emdTab(outI,inI)=EMD2;
        end
    end
    save(['Result/RD/LDM/2/',num2str(i),'/EMD1.mat'],'EMD1');
    save(['Result/RD/LDM/2/',num2str(i),'/T1.mat'],'T1');
    save(['Result/RD/LDM/2/',num2str(i),'/emdTab.mat'],'emdTab');
    save(['Result/RD/LDM/2/',num2str(i),'/timTab.mat'],'timTab');
end