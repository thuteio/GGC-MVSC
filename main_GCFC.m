clc;
clear;
warning('off');
addpath('CVI',"Pre-treatment","FGCG",genpath('GCFC'),"MSA")

%%%%%%%%%%%%%  Data  %%%%%%%%%%%%%%%%%%%
filename = "3Sources";

X=importdata("data/"+filename+".mat");

data=X.data;
label=X.label;

N=size(data{1},1);
V=length(data);
c=length(unique(label));

%%%%%%%%%%%%%  Graph  %%%%%%%%%%%%%%%%%%
switch filename
    case "3Sources"
        p=0.5;  lambda=1e5; theta=1e-4;
    case "100leaves"
        p=0.8; lambda=1e3;theta=1e-2;
    case "MSRC"
        p=1; lambda=1e-1;theta=1e-2;
end

%X.G=[];%%%
%X.GGC=[];%%%

if isempty(X.G)
    data_G=Generate_graph(data,1);
else
    data_G=X.G;
end
%%%%%%%%%%%%%  FGCG  %%%%%%%%%%%%%%%%%%%

if isempty(X.GGC)
    [data_GGC,trans,GGC] = FGCG(data_G,c,p);
    M = MSA(GGC);  
else
    M = X.M;
    trans=X.trans;
    data_GGC=X.GGC;
end

%%%%%%%%%%%%%  Clustering  %%%%%%%%%%%%%%%%%%%

[idx,gamma]=GCFC(data_GGC,c,M,lambda,trans,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ACC,NMI,ARI,purity]=CVI(idx,label);
fprintf(" ACC=%.4f\n NMI=%.4f\n ARI=%.4f\n Purity=%.4f\n=============\n\",ACC,NMI,ARI,purity);