function disp1 = R_IC(kk,ff,disp0)
sdof = length(ff);
sdof0 = length(disp0);
sdof1 = sdof;
%{
if sdof0<sdof
%     disp0(sdof) = 0;
    disp0(sdof0+1:sdof) = kk(sdof0+1:sdof,sdof0+1:sdof)\(ff(sdof0+1:sdof)-kk(sdof0+1:sdof,1:sdof0)*disp0);
elseif sdof0>sdof
    sdof = sdof0;
    kk(sdof,sdof) = 0;
    ff(sdof) = 0;
    for i = sdof1+1:sdof
        kk(i,i) = 1;
    end
end
%}
d = ff-kk*disp0;
vec = abs(d);
% [vec1,ord] = sort(vec,'descend');
index1 = zeros(sdof,1);
n0 = 0;
ep = 1*10^(-8);
for i = 1:sdof
    if vec(i)>ep
        rk = kk(i,:);
        for j = 1:sdof
            if rk(j)~=0&&index1(j)==0
                n0 = n0+1;
                index1(j) = j;
            end
        end
    end
%     if n0>=sdof/10
%         break
%     end
end
n0
index = zeros(n0,1);
n1 = 0;
for i = 1:sdof
    if index1(i)~=0
        n1 = n1+1;
        index(n1) = index1(i);
    end
end

disp1 = disp0;
rb = sparse(sdof,n0);
% rb(:,1) = disp0;
for i = 1:n0
    rb(index(i),i) = 1;
%     rb(index(i-1),1) = 0;
%     disp1(index(i)) = 0;
end

kkr = rb'*kk*rb;
% d = ff-kk*disp1;
ffr = rb'*d;
z = kkr\ffr;
% z1 = z(1)
% disp1 = rb*z;
disp1 = disp1+rb*z;
