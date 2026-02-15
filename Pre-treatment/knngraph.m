function K = knngraph(data,k,string)
[N,~]=size(data);
[nid,dis]=knnsearch(data,data,'K',k+1);
nid=nid(:,2:end);
dis=dis(:,2:end).^2;
gamma=mean(dis,'all');
K=eye([N N]);
for i=1:N
    %dis=pdist2(data(i,:),data(nid(i,:),:),"squaredeuclidean");
    K(i,nid(i,:))=exp(-(dis(i,:)/gamma));
end
if string=="D"
    K=(K+K')/2;
else 
    K=min(K,K');
end

end