function x = GSiter(L,U,x0,ff,nst)

% L = tril(kk,0);
% U = -triu(kk,1);
b = L\ff;

for i = 1:nst
    x1 = U*x0;
    x2 = L\x1;
    x = x2+b;
%     dx = ff-kk*x;
%     i
%     delt = sqrt(dx'*dx)/sqrt(ff'*ff);
%     x0 = smooth(x);
    x0 = x;
end
    