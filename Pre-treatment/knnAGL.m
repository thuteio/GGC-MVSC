function As = knnAGL(data,k)
[N,~]=size(data);
[nid,dis]=knnsearch(data,data,'K',k+2);
dis=dis.^2;
diss=dis(:,end)-dis(:,2:(end-1));
diss=diss./sum(diss,2);
diss(isnan(diss))=0;
As=zeros([N N]);
for i=1:N
    As(i,nid(i,2:end-1))=diss(i,:);
end
As=(As+As')/2;
end