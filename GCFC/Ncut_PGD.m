function [idx,U,alpha,err]=Ncut_PGD(L,D,U,c,maxiter,lambda)
N=size(L,1);
D=diag(D);
gradu=@(u)(2/(u'*D*u))*(L*u-((u'*L*u)/(u'*D*u))*D*u);
Lu=@(u)(u'*L*u)/(u'*D*u);
alpha=1;
gold=(sqrt(5)-1)/2;
for iter=1:maxiter
    Lold=0;
    grad=zeros(size(U));
    for j=1:c
        Lold=Lold+Lu(U(:,j));
        grad(:,j)=gradu(U(:,j));
    end

    while alpha>1e-5
        Unew=U-grad*alpha;
        for i=1:N
            Unew(i,:)=rntproj(Unew(i,:),1e-6,lambda);
        end
        Lnew=0;
        for j=1:c
            Lnew=Lnew+Lu(Unew(:,j));
        end

        if Lnew<Lold
            break
        else
            alpha=alpha*gold;
        end
    end

    err=abs(Lnew-Lold);
    U=Unew;
    if err<1e-4 ||alpha<1e-5
        break
    end

end
[~,idx]=max(U,[],2);
end
