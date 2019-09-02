%one mechanical part1 - the half of Bearing seat 

Xi = [0,0,1,1,2,2,3,3,4,4];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,1];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);


B=cell(8,2,2);
s = 1/sqrt(2);h=2; rr=1;r=0.55;  jj=2*(1-s)/(3-s);

B{1,1,1}=[0   r  0   1];
B{1,1,2}=[0   r  h   1];
B{1,2,1}=[0   rr  0   1];
B{1,2,2}=[0   rr  h   1];


B{2,1,1}=[-r   r  0  s];
B{2,1,2}=[-r   r  h  s];
B{2,2,1}=[-rr   rr   0   s];
B{2,2,2}=[-rr   rr  h    s];

B{3,1,1}=[-r   0  0   1];
B{3,1,2}=[-r   0  h   1];
B{3,2,1}=[-rr   0  0    1];
B{3,2,2}=[-rr   0  h    1];


B{4,1,1}=[-r   -r  0   s];
B{4,1,2}=[-r   -r  h   s];
B{4,2,1}=[-rr   -rr  0  s];
B{4,2,2}=[-rr   -rr  h  s];

B{5,1,1}=[0   -r  0   1];
B{5,1,2}=[0   -r  h   1];
B{5,2,1}=[-rr*s-0.02   -rr*s-0.14  0   1];
B{5,2,2}=[-rr*s-0.02   -rr*s-0.14  h   1];

B{6,1,1}=[0   -rr  0   1];
B{6,1,2}=[0   -rr  h   1];
B{6,2,1}=[-rr+jj   -rr  0   1];
B{6,2,2}=[-rr+jj   -rr  h   1];


B{7,1,1}=[0   -2.5*rr  0   1];
B{7,1,2}=[0   -2.5*rr  h   1];
B{7,2,1}=[-rr  -2.5*rr  0   1];
B{7,2,2}=[-rr  -2.5*rr  h   1];

B{8,1,1}=[0   -2.8*rr  0   1];
B{8,1,2}=[0   -2.8*rr  h   1];
B{8,2,1}=[-rr  -2.8*rr  0   1];
B{8,2,2}=[-rr  -2.8*rr  h   1];

B{9,1,1}=[0   -3*rr  0   1];
B{9,1,2}=[0   -3*rr  h   1];
B{9,2,1}=[-rr  -3*rr  0   1];
B{9,2,2}=[-rr  -3*rr  h   1];


n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
