function disp0 = restatic_CA(ndof,K,F)
'start restatic'
t1 = cputime;

load initial

ndof1 = ndof;

kk1=K;

if ndof0<ndof%���ӽڵ�
    K0(ndof,ndof) = 0;
    disp0(ndof) = 0;
elseif ndof0>ndof%ɾ���ڵ�
    ndof = ndof0;
    K(ndof,ndof) = 0;    %��0�������
    F(ndof) = 0;
    for i = ndof1+1:ndof
        K(i,i) = 1;     %�������Խ��ߴ���1????Ϊɶ
    end
end
deltk = K-K0;


alfa = 0.001;
if ndof0<ndof%���ӽڵ�
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
for i = 2:nb              %%CA�����������rB.
    vec1 = deltk*rb(:,i-1);
    vec2 = inv(K0)*vec1;
    rb(:,i) = -vec2;
end
kkr = rb'*K*rb;
clear kk
ffr = rb'*F;
z = inv(kkr)*ffr;    
disp0 = rb*z;        %%CA����λ��
if ndof1<ndof0  %%ɾ���ڵ�
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