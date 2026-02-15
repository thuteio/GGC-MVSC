function [idx,Fm,gamma]=GCFC(data,c,M,lambda,trans,theta)
[~,~,V]=size(data);
len=size(M,3);

% gamma=rand([V 1]);
% gamma=gamma./(sum(gamma));
gamma=(1/V)*ones([V 1]);
gamma_tensor = reshape(gamma, 1, 1, V);
Xbar = sum(data .* gamma_tensor, 3);
[Fm,snh] = finchpp(Xbar, c);
if snh==-1
    idx=-1;
    gamma=-1;
    return;
end


% [~,idx]=max(full(idx),[],2);
% Fm = full(ind2vec(idx', c))';

eta=zeros([len 1]);
f=zeros([V 1]);
D=sum(data,2);
L=zeros(size(data));
for v=1:V
    L(:,:,v)=diag(D(:,:,v))-data(:,:,v);
end

for iter=1:1e2
    %updata omega
    for i=1:len
        eta(i)=gamma'*M(:,:,i)*gamma;
    end
    omega=1./eta;
    omega=omega./sum(omega);

    %updata gamma
    omega_tensor = reshape(omega.^2, 1, 1, len);
    H=lambda*sum(M .* omega_tensor, 3);
    for v=1:V
        f(v)=trace(Fm'*L(:,:,v)*Fm)/trace(Fm'*diag(D(:,:,v))*Fm);
    end
    gammanew=QPAS(H, f);

    %updata Y
    gamma_tensor = reshape(gammanew, 1, 1, V);
    lbar = sum(L .* gamma_tensor, 3);
    dbar = sum(D .* gamma_tensor, 3);
    [~,Fm,alpha,errl]=Ncut_PGD(lbar,dbar,Fm,c,1,theta);

    %Termination condition
    err=sum(abs(gammanew-gamma));
    gamma=gammanew;
    if err<1e-4||alpha<1e-6||errl<1e-6
        break
    end
end
Fin=trans'*Fm;
[~,idx]=max(Fin,[],2);
end