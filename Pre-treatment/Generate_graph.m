function data_G = Generate_graph(data,f)
N=size(data{1},1);
V=length(data);
data_G=zeros([N N V]);

[normsty,sumD] = G_Parament(data,V,f);
k=20;
for v=1:V
    datal=datanorm(data{v},normsty);
    if sumD>=1e4
        k=20;
        G = knngraph(datal,k,"D");
    else
        if f==0
            k=15;
        end
        G = knnAGL(datal,k);
    end
    G(1:N+1:end) = 0;
    data_G(:,:,v)=G;
end
end