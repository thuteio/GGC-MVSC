function [idxset,NBG]=fgcg_main(A,c,NBG)
A=(A+A')/2;
G=graph(A);
N=size(A,1);

Graset=cell([N 1]);
for i=1:N
    Graset{i}.idx=i;
    Graset{i}.Q=0;
end
idx=(1:N)';
degreeset=sum(A,1);
idxset=[];

while 1
    change=0;
    Nset=randperm(N);
    for iii=1:N
        i=Nset(iii);
        neig=neighbors(G,i);
        inneig = neig(idx(neig)==idx(i));
        Qout=0;
        if ~isempty(inneig)
            idxi=Graset{idx(i)}.idx;
            idxl=setdiff(idxi,i);
            Qoutl=Qscore(degreeset(idxl),A(idxl,idxl));
            Qout=Qoutl-Graset{idx(i)}.Q;
        end

        outneig=setdiff(neig,inneig);
        if ~isempty(outneig)
            Qset=zeros([1 length(outneig)]);
            QinQ=zeros([1 length(outneig)]);
            for j=1:length(Qset)
                ji=outneig(j);
                idxj=Graset{idx(ji)}.idx;
                QinQ(j)=Qscore(degreeset([idxj;i]),A([idxj;i],[idxj;i]));
                Qset(j)=Qout+(QinQ(j)-Graset{idx(ji)}.Q);
            end
            [Qmax,maxidx]=max(Qset);
            if Qmax>0
                change=1;
                Graset{idx(i)}.idx=setdiff(Graset{idx(i)}.idx,i);
                if isempty(Graset{idx(i)}.idx)
                    Graset{idx(i)}.Q=0;
                else
                    Graset{idx(i)}.Q=Qoutl;
                end
                idx(i)=idx(outneig(maxidx));%Add the corresponding module
                Graset{idx(i)}.idx=[Graset{idx(i)}.idx;i];
                Graset{idx(i)}.Q=QinQ(maxidx);
            end
        end
    end
    cci=length(unique(idx));
    if change==0||cci<c
        break
    end
    gbi=0;
    for i=1:length(Graset)
        if isempty(Graset{i}.idx)
            continue;
        else
            gbi=gbi+1;
            idxNBG=Graset{i}.idx;
            NBG(idxNBG,idxNBG)=NBG(idxNBG,idxNBG)+1;
        end
    end
    idxset=[idxset idx];
end
end