function M=msa_main(idxset,f)
V=length(idxset);
M=eye([V V]);
for i=1:V
    for j=(i+1):V
        mi = normalized_DTW_similarity(idxset{i},idxset{j}, f);
        M(i,j)=mi;
        M(j,i)=mi;
    end
end
end


function similarity = normalized_DTW_similarity(s1, s2, f)


    m = size(s1,2);
    n = size(s2,2);
    dp = zeros(m, n);
    dp(1,1)=f(s1(:,1),s2(:,1));
    path_length = zeros(m, n);
    path_length(1,1)=1;
    for i=2:m
        dp(i,1)=dp(i-1,1)+f(s1(:,i),s2(:,1));
        path_length(i,1)=path_length(i-1,1)+1;
    end
    for j=2:n
        dp(1,j)=dp(1,j-1)+f(s1(:,1),s2(:,j));
        path_length(1,j)=path_length(1,j-1)+1;
    end

    for i=2:m
        for j=2:n
            can=[dp(i-1,j-1),dp(i,j-1),dp(i-1,j)];
            len=[path_length(i-1,j-1),path_length(i,j-1),path_length(i-1,j)];
            [~,idx]=max(can);
            dp(i,j)=f(s1(:,i),s2(:,j))+can(idx);
            path_length(i,j)=1+len(idx);
        end
    end
    similarity=dp(m,n)/path_length(m,n);

end