function x = AMG111(kk,ff,kk11,P,p1,x0)

sdof1 = length(kk11(:,1));
sdof0 = length(p1(:,1));
% if(sdof0<sdof1)
%     k11 = kk11(1:sdof0,1:sdof0);
%     k12 = kk11(1:sdof0,sdof0+1:sdof1);
%     k21 = kk11(sdof0+1:sdof1,1:sdof0);
%     k22 = kk11(sdof0+1:sdof1,sdof0+1:sdof1);
% %     ik = inv(full(k22));
% %     km = k12*ik;
%     km=k12/k22;
%     km = sparse(km);
%     kk11 = k11-km*k21;
% end
kk11=sparse(kk11);
kk1 = p1'*kk11*p1;

 for i = 1:size(kk1,1)
    if(abs(kk1(i,i))<1e-6)
        kk1(i,i) = 1;
    end
 end
% i1=0:size(kk1,1)-1;
% ib=i1*(size(kk1,1))+i1+1;
% %  i1 = 1:size(kk1,1);
%     if(abs(kk1(ib))<1e-6)
%         kk1(ib) = 1;
%     end
% kk2 = inv(kk1);



% kk2 = p2'*kk1*p2;
% % ik2 = inv(kk1);
%  for i = 1:size(kk2,1)
%     if(abs(kk2(i,i))<1e-6)
%         kk2(i,i) = 1;
%     end
%  end

% ik2 = inv(kk2);
L1 = tril(kk11,0);
U1 = -triu(kk11,1);  
L = tril(kk,0);
U = -triu(kk,1);
clear kk0
n = length(ff);
% x0 = zeros(n,1);
x0 = P*x0;
x = GSiter(L,U,x0,ff,3);  %前光滑
 for i = 1:10
d = ff-kk*x;
% i
deltd = sqrt(d'*d)/sqrt(ff'*ff);
dx = x-x0;
deltx = sqrt(dx'*dx)/sqrt(x'*x);
% if(mod(i,10)==1)
%     i
%     deltd
%     deltx
% end
if(deltx<1e-20)
    deltx
    break;
end


% if(deltd<1e-4)
%     deltd
%     break;
% end
d1 = P'*d;


% dx1 = CA(invkk,kk1,dk,d1);
% dx1 = lanczos(invkk,kk1,dk,d1);
% dx1 = kk1\d1;

% dx1 = invkk*d1;

if(sdof0<sdof1)
    f1 = d1(1:sdof0);
    f2 = d1(sdof0+1:sdof1);
    d1 = f1-km*f2;
end
dx1 = zeros(length(d1),1);
dx1 = GSiter(L1,U1,dx1,d1,3);  %：使用初始解（x0=0）对方程应用光滑算子如式（3.75）  
% dx1 = iteration(kk0,kk1,ik2,p1,p2,d1);
% dx1 = iteration(kk11,kk1,ik2,p1,p2,d1);   %重分析所在
%% 残余向量
d1 = d1-kk11*dx1;   %1

%% 限制残余向量到网格M2
d2 = p1'*d1;      %2    

%% 算法2 第5步，MG计算v(l+1)
% dx2 = kk2\d2;
dx2 = kk1\d2;   %用MG计算v(l+1),如算法2中的第5步
dx = p1*dx2;   
dx1 = dx1+dx;  

dx1 = GSiter(L1,U1,dx1,d1,3);  


if(sdof0<sdof1)
    d2 = f2-k21*dx1;
    dx2 = k22\d2;
    dx1 = [dx1;dx2];
end
       
% dx1 = direct(invkk,kk1,dk,d1);
dx = P*dx1;
x0 = x+dx;

x = GSiter(L,U,x0,ff,3);

 end
%     i
%     deltd
%     deltx
