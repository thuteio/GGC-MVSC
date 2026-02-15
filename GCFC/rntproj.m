function [s] = rntproj(t,th,lambda)
%初始化
a=mean(t);
if isnan(a)
    t=datanorm(t,'range');
    a=mean(t);
end
while 1
ft=sum(max((t+a)/(1+lambda),0))-1;
dt=sum(max(sign((t+a)/(1+lambda)),0));
if dt==0
    dt=1e-4;
end
anew=a-ft/dt;
if abs(anew-a)<=th
    break;
end
a=anew;
end
s=max((t+a)/(1+lambda),0);
end