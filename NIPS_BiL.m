clear
clc
%% Bilinguistics
filefds=dir('data/');
Num=numel(filefds);% the number of language pairs
for i=3:Num
    i-2
    Name=filefds(i).name;% the name of file folder
    load(['data/',Name,'/A.mat']);
    load(['data/',Name,'/B.mat']);
    load(['data/',Name,'/WA.mat']);
    load(['data/',Name,'/WB.mat']);
    if Name=='zh-en'
        IT=load(['data/',Name,'/W']);
        IT=IT';% the initialization of transformation matrix
    else
        IT=randn(size(A,1));
    end
    %% Original EMD
    disp('Original EMD');
    tic
    [EMD1,~,~]=SinkhornInit(A,B,WA,WB,IT);
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
        k1=ceil(size(A,2)/csArr(outI));% the size of core-set of set A
        k2=ceil(size(B,2)/csArr(outI));% the size of core-set of set B
        for inI=1:runNo
            outI,inI
            tic
            [CSA,CSWA]=KCenter(A,WA,k1);
            [CSB,CSWB]=KCenter(B,WB,k2);
            [~,~,TM2]=SinkhornInit(CSA,CSB,CSWA,CSWB,IT);
            [EMD2,~]=Sinkhorn(TM2*copA,copB,copWA,copWB);
            timTab(outI,inI)=toc;
            emdTab(outI,inI)=EMD2;
        end
    end
    save(['Result/BiL/',Name,'/EMD1.mat'],'EMD1');
    save(['Result/BiL/',Name,'/T1.mat'],'T1');
    save(['Result/BiL/',Name,'/emdTab.mat'],'emdTab');
    save(['Result/BiL/',Name,'/timTab.mat'],'timTab');
end