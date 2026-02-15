function [idx,U]=Ncut_CD(A,D,U,c)

N=size(A,1);
%Dll=diag(sum(A,2));
D=diag(D);
uTA = U'*A; 

uTAu = zeros(c, 1);
uTDu = zeros(c, 1);
for j = 1:c
    uTAu(j) = uTA(j,:) * U(:, j); 
    uTDu(j) = U(:, j)' * D * U(:, j); 
end

for iter=1:1e2
    change=false;
    for i=1:N
        psj=zeros([1 c]);
        [~,p]=max(U(i,:),[],2);
        for j=1:c
            temp1=uTAu(j)/uTDu(j);
            if j==p
                temp2=(uTAu(j)-2*uTA(j,i)+A(i,i))/(uTDu(j)-D(i,i));
                psj(j)=temp1-temp2;
            else
                temp2=(uTAu(j)+2*uTA(j,i)+A(i,i))/(uTDu(j)+D(i,i));
                psj(j)=temp2-temp1;
            end
        end
        [~, r] = min(psj);
        if r~=p
            change=true;
            %U
            U(i,p)=0;
            U(i,r)=1;

            %uTAu
            uTAu(p)=U(:, p)' * A * U(:, p);
            uTAu(r)=U(:, r)' * A * U(:, r);

            %uTA
            uTA(p,:)=uTA(p,:)-A(i,:);
            uTA(r,:)=uTA(r,:)+A(i,:);
            
            %uTDu
            uTDu(p)=U(:, p)' * D * U(:, p);
            uTDu(r)=U(:, r)' * D * U(:, r);
        end
    end
    if ~change
        break;
    end
end
[~,idx]=max(U,[],2);
end