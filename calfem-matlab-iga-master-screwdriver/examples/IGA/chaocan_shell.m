clear all; close all; clc;

%% Load data   (3个自由度的solid单元)
addpath(genpath('C:\Users\Administrator\Desktop\IGA\calfem-matlab-iga-master'));

% Read input file
geometry=geo_read_nurbs( 'geo_beam3D.txt');

%--------------------------------------------------------------------------
% %plot original NURBS volum and control points

%figure(1)
%nrbkntplot (geometry.nurbs);
 %hold on;
%nrbctrlplot(geometry.nurbs);

%--------------------------------------------------------------------------
degree     = [2 1 1];     
regularity = [1 1 1];     
nquad      = [4 4 4];     
nsub       = [10 0 1];     
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);
KV.Xi =nurbs.knots{1}; KV.Eta = nurbs.knots{3}; KV.Zeta = nurbs.knots{2}; clear Xi Eta Zeta
 deg.p = degree(1); deg.q=degree(3); deg.r = degree(2); clear p q r
B=cell(nurbs.number(1),nurbs.number(3));
for i=1:nurbs.number(1)
    for j=1:nurbs.number(3)
B{i,j}=nurbs.coefs(:,i,2,j)';
    end
end
[n,m]=size(B)
l=size(B,3)
%{
 su_;

% run ./../NURBS/Surface/ex4
 %run ./../NURBS/solid/data_solid7

%判断二维或三维


 deg.p = p; deg.q=q; deg.r = r; clear p q r


% Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
[ X_ ] = getSubDivKVValues( Xi, 3);
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta, 1);
[B,Eta] = RefineKnotSolidEta( B,deg.q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta, 1);
% [B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );

% Number of control points after refinement
n = size(B,1);
m = size(B,2);
l = size(B,3);

%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta
%}

% Build connectivity arrays
%nel = (n-deg.p) * (m-deg.q) * (l-deg.r);       % number of elements-三维
nel = (n-deg.p) * (m-deg.q);

nnp = n*m; % number of global basis functions
% nen = (deg.p+1)*(deg.q+1)*(deg.r+1);           % number of local basis functions
nen = (deg.p+1)*(deg.q+1);

ndof = nnp*5; % number of global degrees of freedom
ldof = nen*5; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg,n,m,l ); % = INn (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN))*5,5,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(5*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed)
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),5*nen,1);
end

%% Material parameters:
E_Y = 1e4;   ti=1;
nu_P = 0.3;
denisty=785;
damp = 0.1;
D=hooke_strain(4,E_Y,nu_P);

%% Gauss-Legendre quadrature points:
[ gp_x,w_x ] = getGP( deg.p );
[ gp_y,w_y ] = getGP( deg.q );
[ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation

% Element loop
K = zeros(ndof); % Needs to be changed to sparse for large problems!!
M = zeros(ndof);
C = zeros(ndof);
F = zeros(ndof,1);
for e = 1 : nel
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1)
    nj = INN(IEN(1,e),2)
%    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) 
        continue
    end
    
    Ke = zeros(nen*5);
    Fe = zeros(nen*5,1);
    Me = zeros(nen*5);
    Ce = zeros(nen*5);jishu=0;
    for i = 1 : NQUADx % Loop trough Gauss points高斯积分
        for j = 1 : NQUADy
          for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                jishu=jishu+1
                
                
                %%---------------超参壳元（mindlin壳）
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ Bm,J,R ] = Shape_function_modify( GP,e,deg,B,KV,INN,IEN);
                J
                % Combine quadrature weights with det(J)
                Jmod = abs(det(J))*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                Ke_ = Bm'*D*Bm*Jmod;
%                 [ Ke_ ] = Build_K_Local( dR_dx,Jmod,D,nen );   %solid 
               % [ Ke_ ] = Build_K_Local_shell( dR_dx,Jmod,D,nen );   %超参shell 
              
                
                %{
                %%------***-------kirchhoff-love 壳---------
                
                [ Bm,Bb,J,R,T ] = Shape_function_modify( GP,e,deg,B,KV,INN,IEN); 
                Jmod = abs(J)*w_x(i)*w_y(j)*w_z(k);  %
                Km = ti*Bm'*D*Bm*Jmod;
               Kb = ti^3/12* Bb'*D*Bb*Jmod; 
               Kee = Km+Kb;
                  
%                 Kee =unit_assemble_stiff(Km,Kb,nen);
                %} 
             
                Ke = Ke + Ke_;
%                  Ke = Ke + Kee;
                %{
                % Build Me
                 Ni=eye(5);
                 for s=1:length(R)
                     shamtx(:,(s-1)*5+1:s*5)=R(s)*Ni;
                 end
                 Me_=denisty*shamtx'*shamtx*Jmod;
                 Ce_=damp*shamtx'*shamtx*Jmod;
                 Me = Me + Me_;
                 Ce = Me + Me_;
                %}
          end
        end
    end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
    M(idx,idx) = M(idx,idx) + Me;
    C(idx,idx) = C(idx,idx) + Ce;
    F(idx) = F(idx)+Fe;
end


% Apply load at two corner nodes
loc_z=1;
loc_y=8;
constNod = [];
for i = 1 : numel(B)
    if B{i}(3) == loc_z
        if B{i}(1) == 48
            if B{i}(2) == loc_y
                constNod=[constNod i];
            end
        end
        
    end
end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:
F(ID(3,constNod)) =1;

%% Boundary condiitons

% Find controlpoints with x = 0 to constrain
loc_x=0;
constNod = [];
for i = 1 : numel(B)
    if B{i}(1) == loc_x
        constNod=[constNod i];
    end
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];

for i = 1:ndof
    rk = K(i,:);
    
    k = find(abs(rk)>0); 
    if length(k)==0
        K(i,i) = 1.0e-10;      %将距阵中全为0的行元素极小化
    end
  
end

% loc_f=4;
% constf = [];
% for i = 1 : numel(B)
%     if B{i}(1) == loc_f
%         constf=[constf i];
%     end
% end
% bf1=reshape(ID(:,constf),numel(ID(:,constf)),1);
% bc=[bc1,zeros(length(bc1),1)];



%% Solve system
[a,K]=solveq(K,F,bc);

%{
%% 时域初始计算 Newmark法

%给定初始值

deltt=0.01;alf=0.25;delt=0.5;
c0=1/(alf*deltt^2);c1=delt/(alf*deltt);c2=1/(alf*deltt); c3=1/(2*alf)-1;
c4=delt/alf-1; c5=deltt/2*(delt/alf-2); c6=deltt*(1-delt);  c7=delt*deltt;

%形成有效的质量矩阵
Md=K+c0*M+c1*C;
U=chol(Md);

dt= zeros(ndof,1); vt = zeros(ndof,1); at = zeros(ndof,1); att = zeros(ndof,1); t=0; Q=F;
for n=1:1/deltt
   t=t+deltt;
   %claculate the payload 计算有效载荷
   Qe=Q+M*(c0*dt+c2*vt++c3*at)+C*(c1*dt+c4*vt+c5*at);
    %claculate the displacement  U*U'*dtt= Qe
    mid=U\Qe;
    dtt=U'\mid;
     %claculate the velocity and acceleration
    att=c0*(dtt-dt)-c2*vt-c3*att;
    for j=1:length(bf1)
    vt(j,1)=5000;   %常速度载荷
    end
    %赋给下次迭代 
    dt=dtt;
    at=att;
    
end
%}
% Rewrite to cell structure. (For plotting)
%{
u = cell(size(B));
comb = cell(size(B));
for i = 1 : size(ID,2)
    u{i} = [dt(ID(:,i)); 0];
    comb{i} = B{i} + u{i};
end
%}

u = cell(size(B));
comb = cell(size(B));
for i = 1 : size(ID,2)
    u{i} = [a(ID(1:3,i))',0];
    comb{i} = B{i} + u{i};
end
for i=1:length(B)
Sx(i,1)=B{i,1}(1);SX(i,1)=comb{i,1}(1);
Sy(i,1)=B{i,1}(2);SY(i,1)=comb{i,1}(2);
Sz(i,1)=B{i,1}(3);SZ(i,1)=comb{i,1}(3);
Sx(i,2)=B{i,2}(1);SX(i,2)=comb{i,2}(1);
Sy(i,2)=B{i,2}(2);SY(i,2)=comb{i,2}(2);
Sz(i,2)=B{i,2}(3);SZ(i,2)=comb{i,2}(3);
end
figure(1)
surf(Sx,Sy,Sz)
shading interp
colormap summer
hidden off
axis equal
title('NURBS surface')
xlabel('x')
ylabel('y')
zlabel('z')
hold on


figure(2)
surf(SX,SY,SZ)
shading interp
colormap summer
hidden off
axis equal
title('NURBS surface')
xlabel('x')
ylabel('y')
zlabel('z')
hold on
%% plot
% plotNurbsSolidElementSimple( KV, B )
% title('Initial geometry')
% plotNurbsSolidElementSimple( KV, comb )
% title('Displacements')