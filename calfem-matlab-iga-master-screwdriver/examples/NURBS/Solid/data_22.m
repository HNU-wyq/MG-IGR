%one mechanical part1 -Bearing seat 

Xi =[0,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
% Zeta = [0,1,2,3,4,5,6,7,8];
Zeta = [0,0,1,1,2,3,4,4,5,5];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);



B=cell(2,2,9);
s = 1/sqrt(2); r3=0.3; rr=1; rx=-1.75;rz=0.5;
s1=(rx+r3)/(rx+r3*s);s2=(rx-r3*s)/(rx-r3);s3=(rz-r3)/(rz-r3*s);s4=(rz+r3*s)/(rz+r3);

B{2,1,1}=[rx+r3*s   -2.5*rr  rz-r3*s   1];
B{2,1,2}=[rx   -2.5*rr   rz-r3   1];
B{2,1,3}=[rx-r3*s   -2.5*rr   rz-r3*s   1];
B{2,1,4}=[rx-r3   -2.5*rr   rz   1];
B{2,1,5}=[rx-r3*s   -2.5*rr   rz+r3*s   1];
B{2,1,6}=[-1.75   -2.5*rr   rz+r3   1];
B{2,1,7}=[rx+r3*s   -2.5*rr   rz+r3*s   1];
B{2,1,8}=[rx+r3   -2.5*rr   rz   1];
B{2,1,9}=[rx+r3*s   -2.5*rr  rz-r3*s   1];

B{2,2,1}=[rx+r3*s   -3*rr  rz-r3*s   1];
B{2,2,2}=[rx   -3*rr   rz-r3   1];
B{2,2,3}=[rx-r3*s   -3*rr   rz-r3*s   1];
B{2,2,4}=[rx-r3   -3*rr   rz   1];
B{2,2,5}=[rx-r3*s   -3*rr   rz+r3*s   1];
B{2,2,6}=[-1.75   -3*rr   rz+r3   1];
B{2,2,7}=[rx+r3*s   -3*rr   rz+r3*s   1];
B{2,2,8}=[rx+r3   -3*rr   rz   1];
B{2,2,9}=[rx+r3*s   -3*rr  rz-r3*s   1];

B{1,1,1}=[-1   -2.5*rr   0   1];
B{1,1,2}=[-1.75   -2.5*rr   0   1];
B{1,1,3}=[-2.5   -2.5*rr   0   1];
B{1,1,4}=[-2.5   -2.5*rr   0.5   1];
B{1,1,5}=[-2.5   -2.5*rr   1   1];
B{1,1,6}=[-1.75   -2.5*rr   1   1];
B{1,1,7}=[-1   -2.5*rr   1   1];
B{1,1,8}=[-1   -2.5*rr   0.5   1];
B{1,1,9}=[-1   -2.5*rr   0   1];

B{1,2,1}=[-1   -3*rr   0   1];
B{1,2,2}=[-1.75   -3*rr   0   1];
B{1,2,3}=[-2.5   -3*rr   0   1];
B{1,2,4}=[-2.5   -3*rr   0.5   1];
B{1,2,5}=[-2.5   -3*rr   1   1];
B{1,2,6}=[-1.75   -3*rr   1   1];
B{1,2,7}=[-1   -3*rr   1   1];
B{1,2,8}=[-1   -3*rr   0.5   1];
B{1,2,9}=[-1   -3*rr   0   1];

n = size(B,1);m = size(B,2);l = size(B,3);dim=4;
Bm=cell2mat(B);
trans =[
    1.0000   0.0000         0         0
    0.0000    1.0000         0         0
    0         0    -1.0000         0
    0         0         0    1.0000];
Bm2=permute(Bm,[2,1,3]);
ff=reshape(Bm2,dim,n*m*l);ff1=num2cell(ff',[ 2 3]);

Bm1 =num2cell( (trans*ff)',[ 2 3]);

solid12.coefs=reshape(Bm1,[2,2,9]);
% Order of basis
p = k_-n-1;q = l_-m-1;r = m_ -l -1;

solid11=nrbmak(B,{Xi Eta Zeta});solid12=solid11;
solid12.coefs=reshape(Bm1,[2,2,9]);
% BB=cat(3,solid11.coefs{:});
% BB(:,3,:)=BB(:,3,:)+1;
% solild12.coefs=cellfun(@ ,solid12.coefs);

numpatches=2;
patches(1)=solid11;patches(2)=solid12;
%% order the number of node and element of patches 
for ip=1:numpatches
    B=patches(ip).coefs;
% B=solid11.coefs;
   BJ=cell2mat(B);
   Xi = patches(ip).knots{1};
   Eta = patches(ip).knots{2};
   Zeta= patches(ip).knots{3};
% n1=size(B,1);m1=size(B,2);l1=size(B,3);
% [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );

%    B=B(1:20,:,:);
% Number of control points after refinement
n = size(B,1);
m = size(B,2);
l = size(B,3);

%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) * (l-deg.r); % number of elements
nnp = n*m*l; % number of global basis functions

nen = (deg.p+1)*(deg.q+1)*(deg.r+1); % number of local basis functions
ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
 for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
 end
 
loc_z=0;
commonNod = [];
 for i = 1 : numel(B)
        if B{i}(3) == loc_z
                commonNod=[commonNod i];
        end
 end
 
 Inn(ip)=INN;Ien(IP)=IEN; IID(ip)=ID;Lm(ip)=LM;
 ComNod(ip)=commomNod;
end

%% reorder the number of the node of patches except the first patch

