clear
clc
%% PPI
filefds=dir('data/PPI/emb2/*.mat');
Num=numel(filefds);% the number of networks
data=cell(1,Num);
fileNms=cell(1,Num);
for i=1:Num
    tp1=load(['data/PPI/emb2/',filefds(i).name]);
    tp2=fieldnames(tp1);
    fileNms{i}=tp2{1};
    data{i}=tp1.(tp2{1});
end
clear tp1
csArr=[50,40,30,20,10];% the array of core-set size
csSz=numel(csArr);% the number of core-sets
runNo=20;
for i=1:Num-1
    for j=i+1:Num
        i,j
        A=data{i};
        B=data{j};
        dim=size(A,1);
        n1=size(A,2);% the number of points in set A
        n2=size(B,2);% the number of points in set B
        WA=rand(1,n1);%ones(1,n1);
        WB=rand(1,n2);%ones(1,n2);
        IT=randn(dim);
        %% Original data
        tic
        [EMD1,~,~]=SinkhornInit(A,B,WA,WB,IT);
        T1=toc;
        %% Core-sets
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
            end
        end
        save(['Result/PPI/',...
            fileNms{i},'-',fileNms{j},'/EMD1.mat'],'EMD1');
        save(['Result/PPI/',...
            fileNms{i},'-',fileNms{j},'/T1.mat'],'T1');
        save(['Result/PPI/',...
            fileNms{i},'-',fileNms{j},'/emdTab.mat'],'emdTab');
        save(['Result/PPI/',...
            fileNms{i},'-',fileNms{j},'/timTab.mat'],'timTab');
    end
end