
%%����һ��

Xi = [0,0,1,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,0,1,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,0,1,2,2];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);



B=cell(3,3,4);
s=cos(0.6); h=3;  %
B{1,1,1}=[0     0  0  1];
B{1,1,2}=[-0.7  4.5  0  1];
B{1,1,3}=[0     4.7  0  1];
B{1,1,4}=[0.8   8  0  1];

B{1,2,1}=[1.9   1  0  s];
B{1,2,2}=[1.9   5.5  0  s];
B{1,2,3}=[1.9   5.7  0  s];
B{1,2,4}=[1.9   9  0  s];

B{1,3,1}=[3.8   0  0  1];
B{1,3,2}=[4.5   4.5  0  1];
B{1,3,3}=[3.8   4.7  0  1];
B{1,3,4}=[3     8  0  1];

B{2,1,1}=[0     0    h/2   1];
B{2,1,2}=[-0.7  4.5  h/2  1];
B{2,1,3}=[0     4.7  h/2  1];
B{2,1,4}=[0.8   8  h/2  1];

B{2,2,1}=[1.9   1   h/2  s];
B{2,2,2}=[1.9   5.5  h/2  s];
B{2,2,3}=[1.9   5.7  h/2  s];
B{2,2,4}=[1.9   9  h/2  s];

B{2,3,1}=[3.8   0    h/2  1];
B{2,3,2}=[4.5   4.5  h/2  1];
B{2,3,3}=[3.8   4.7  h/2  1];
B{2,3,4}=[3     8   h/2  1];

B{3,1,1}=[0     0    h  1];
B{3,1,2}=[-0.7  4.5  h  1];
B{3,1,3}=[0     4.7  h  1];
B{3,1,4}=[0.8   8   h  1];

B{3,2,1}=[1.9   1    h  s];
B{3,2,2}=[1.9   5.5  h  s];
B{3,2,3}=[1.9   5.7  h  s];
B{3,2,4}=[1.9   9    h  s];

B{3,3,1}=[3.8   0    h  1];
B{3,3,2}=[4.5   4.5  h  1];
B{3,3,3}=[3.8   4.7  h  1];
B{3,3,4}=[3     8    h  1];

n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
