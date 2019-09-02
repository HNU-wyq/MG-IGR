function x = AMG122(kk,ff,kk11,P,x0)

sdof0 = length(kk11(:,1));
sdof1 = size(P,1);
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
% kk1 = p1'*kk11*p1;
% kk2 = p2'*kk1*p2;

% ik2 = inv(kk11);

    
L = tril(kk,0);
U = -triu(kk,1);
clear kk0
% n = length(ff);
% x0 = zeros(n,1);
x0 = P*x0;
x = GSiter(L,U,x0,ff,3);  %前光滑
for i = 1:70
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
% if(deltd<1e-20)
%        deltd;
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
   
% % dx1 = iteration(kk0,kk1,ik2,p1,p2,d1);
%  dx1 = iteration(kk11,kk1,ik2,p1,p2,d1);   %重分析所在
dx1 = kk11\d1;
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
