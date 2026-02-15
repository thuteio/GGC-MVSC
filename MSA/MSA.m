function M = MSA(GGC)
V=length(GGC);
M=zeros([V V 3]);
M(:,:,1)=msa_main(GGC,@(x,y)ARI(x,y));
M(:,:,2)=msa_main(GGC,@(x,y)NMI(x,y));
M(:,:,3)=msa_main(GGC,@(x,y)purity(x,y));
end