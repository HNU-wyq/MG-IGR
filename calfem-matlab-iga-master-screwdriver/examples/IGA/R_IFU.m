function [a1,dx,d1,d2,d,Fr,sd] = R_IFU(K,F,L,K0,disp0)
tic

load initial
ndof = length(F);
ndof0 = length(disp0);
sdof1 = ndof;
if ndof0<ndof
    for i = ndof0+1:ndof
        K0(i,i) = 1;
        L(i,i) = 1;
        disp0(i) = F(i);
    end
elseif ndof0>ndof
    ndof = ndof0;
    K(ndof,ndof) = 0;
    F(ndof) = 0;
    for i = sdof1+1:ndof
        K(i,i) = 1;
    end
end

%record unbalanced DOFs
% d1 = abs(diag(K)-diag(K0));
 d1 = sum(abs(K-K0),2);
d = F-K*disp0;
d2 = d1+abs(d);
sd = [];
%{
for i = 1:sdof
    if abs(d1(i))>1e-6
        sd = [sd,i];
    end
end
nd = length(sd)
%}
for i = 1:ndof
    if abs(d2(i))>1e-7
        sd = [sd,i];
    end
end
nd = length(sd)

%calculate the right-hand vecrots or extra constrains
Rb = zeros(ndof,nd);
Rb1 = zeros(length(sd),ndof);
Fr = zeros(length(sd),1);
for i = 1:length(sd)
    
    Rb(:,i) = -K(:,sd(i));
    Rb1(i,:) = K(sd(i),:);
    Fr(i) = d(sd(i));
    d(sd(i)) = 0;
end
for i = 1:nd
    Rb(sd(i),:) = 0;
    Rb(sd(i),i) = 1;
end

%apply extra constrains
V = sparse(ndof,nd);
tl1 = cputime;
L = full(L);
for i = nd:-1:1

    l = L(:,sd(i));
    l(sd(i)) = 0;
    V(:,i) = l;
    
%     ri = find(abs(L(:,sd(i)))>0);
%     cj = sd(i)*ones(length(ri),1);
%     su = sub2ind(size(L),ri,cj);
%     L(su) = 0;
%     cj = find(abs(L(sd(i),:))>0);
%     ri = sd(i)*ones(1,length(cj));
%     su = sub2ind(size(L),ri,cj);
%     L(su) = 0;
%     ri = sd(i);
%     cj = sd(i);
%     su = sub2ind(size(L),ri,cj);
%     L(su) = 1;

    L(:,sd(i)) = 0;
    L(sd(i),:) = 0;
    L(sd(i),sd(i)) = 1;
end
tl2 = cputime;
dtl = tl2-tl1;

%==============SMW formula================

B = [V,Rb,d];
 %parfor i = 1:length(B(1,:))
%    B1(:,i) = L\B(:,i);
%    U = L';
%     B(:,i) = U\B1(:,i);
% end

B1 = L\B;
L = L';
B = L\B1;
T = B(:,1:nd);
X0 = B(:,nd+1:2*nd+1);

clear B
    
K = eye(nd);
K = K+V'*T;
V1 = V'*X0;
V2 = K\V1;
clear V1
dX = T*V2;
clear V2
X = X0-dX;
B = X(:,1:nd);
x0 = X(:,nd+1);
clear X0 dX L V K

%{
B = [V,Rb];
B1 = L\B;
L = L';
B = L\B1;
T = B(:,1:nd);
X0 = B(:,nd+1:2*nd);
clear B
K = eye(nd);
K = K+V'*T;
V1 = V'*X0;
V2 = K\V1;
clear V1
dX = T*V2;
clear V2
B = X0-dX;
Kr=Rb1*B;
z = Kr\Fr;
dx = B*z;
a1 = disp0+dx;
if sdof1<ndof0  %%删除节点
    a1 = a1(1:sdof1); 
    ndof = sdof1;
end
%}
%=====================================

Fr = Fr-Rb1*x0;
cond(Fr)
Kr = Rb1*B;
z = Kr\Fr;
dx = B*z+x0;
 
a1 = disp0+dx;
if sdof1<ndof0  %%删除节点
    a1 = a1(1:sdof1); 
    ndof = sdof1;
end

toc
