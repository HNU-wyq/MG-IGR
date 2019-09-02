function disp0 = restatic_CA(ndof,K,F)
'start restatic'
t1 = cputime;

load initial

ndof1 = ndof;

kk1=K;

if ndof0<ndof%增加节点
    K0(ndof,ndof) = 0;
    disp0(ndof) = 0;
elseif ndof0>ndof%删除节点
    ndof = ndof0;
    K(ndof,ndof) = 0;    %补0扩充距阵，
    F(ndof) = 0;
    for i = ndof1+1:ndof
        K(i,i) = 1;     %扩充矩阵对角线处补1????为啥
    end
end
deltk = K-K0;


alfa = 0.001;
if ndof0<ndof%增加节点
    deltk1 = deltk(ndof0+1:ndof,ndof0+1:ndof);
    deltk2 = deltk(ndof0+1:ndof,1:ndof0);
    invdk = full(inv(deltk1));
    disp0(ndof0+1:ndof) = invdk*(F(ndof0+1:ndof)-deltk2*disp0(1:ndof0));
    deltk(ndof0+1:ndof,ndof0+1:ndof) = (1-alfa)*deltk1;
%     deltk1 = alfa*deltk1;
%     invdk = full(inv(deltk1));
    invdk = invdk/alfa;
    invkk(ndof,ndof) = 0;
    invkk(ndof0+1:ndof,ndof0+1:ndof) = invdk;
end
 nb = 10;
rb = zeros(ndof,nb);
rb(:,1) = disp0;
for i = 2:nb              %%CA法构造基向量rB.
    vec1 = deltk*rb(:,i-1);
    vec2 = inv(K0)*vec1;
    rb(:,i) = -vec2;
end
kkr = rb'*K*rb;
clear kk
ffr = rb'*F;
z = inv(kkr)*ffr;    
disp0 = rb*z;        %%CA法求位移
if ndof1<ndof0  %%删除节点
    disp0 = disp0(1:ndof1); 
    ndof = ndof1;
end
save reana disp0 
t2 = cputime;
dt(4) = t2-t1;
'restatic complete'
save reanatime dt
i = [1:ndof:ndof];
dX = disp0(i);
dY = disp0(i+1);
save dis dX dY