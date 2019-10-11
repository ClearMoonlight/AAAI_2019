clear
clc
%% Ground truth
load('data/es-en/en.mat');
load('data/es-en/es.mat');
load('data/es-en/es_en.mat');
%% Test
load('Result/es-en/FM1.mat');
load('Result/es-en/FM.mat');
%% Preprocessing
[~,lef]=setdiff(es_en(:,2),string(''));
es_en=es_en(sort(lef),:);
%% First
count1=0;
ct1=0;
num1=4774;
num2=6637;
k=10;
for i=1:num1
    pos=find(es(i)==es_en(:,1));
    sz=numel(pos);
    if sz>0
        tolStr='';
        temp=FM1(:,i);
        [~,seq]=sort(temp,'descend');
        seq=seq(1:k);
        for j=1:sz
            tolStr=[tolStr,char(es_en(pos(j),2))];
        end
        for j=1:num2
            if strfind(tolStr,en(j))
                break;
            end
        end
        if j<=num2
            count1=count1+1;
        else
            continue;
        end
        for j=1:k
            if strfind(tolStr,en(seq(j)))
                ct1=ct1+1;
                break;
            end
        end
    end
end
f1=ct1/count1
%% Core-set
f1Tab=zeros(5,20);
for i1=1:5
    FMi=FM{i1};
    for j1=1:20
        fmi=FMi(:,:,j1);
        count2=0;
        ct2=0;
        for i=1:num1
            pos=find(es(i)==es_en(:,1));
            sz=numel(pos);
            if sz>0
                tolStr='';
                temp=fmi(:,i);
                [~,seq]=sort(temp,'descend');
                seq=seq(1:k);
                for j=1:sz
                    tolStr=[tolStr,char(es_en(pos(j),2))];
                end
                for j=1:num2
                    if strfind(tolStr,en(j))
                        break;
                    end
                end
                if j<=num2
                    count2=count2+1;
                else
                    continue;
                end
                for j=1:k
                    if strfind(tolStr,en(seq(j)))
                        ct2=ct2+1;
                        break;
                    end
                end
            end
        end
        f1Tab(i1,j1)=ct2/count2;
    end
end