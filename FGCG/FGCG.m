function [data,trans,GGC] = FGCG(data,c,p)
[N,~,V]=size(data);
GGC=cell([V,1]);
NBG=zeros([N,N],"int8");
totp=0;


for v=1:V
    [GGC{v},NBG]=fgcg_main(data(:,:,v),c,NBG);
    totp=totp+size(GGC{v},2);
end


%Space
NBG(1:N+1:end) = 0;
trun=round(totp*p);
NBG=NBG>trun;
NBG_G=graph(NBG);
idxNBG=conncomp(NBG_G);
gnum=length(unique(idxNBG));
trans=full(ind2vec(idxNBG, gnum));
D=sum(trans,2).^(-1/2);
transD=D.*trans;
data=pagemtimes(pagemtimes(transD,data),transD');
for v=1:V
    datal=data(:,:,v);
    Ni=size(datal,1);
    datal(1:1+Ni:end)=0;
    data(:,:,v)=datal;
end

end