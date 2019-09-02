tic
%% Load data
addpath(genpath('C:\Users\Administrator\Desktop\calfem-matlab-iga-master'));
count=1;
%{
%  data_l;
%      run ./../NURBS/Solid/data_solid_chilun
%          run ./../NURBS/solid/data_solid10
          run ./../NURBS/Solid/data_part2
%             run ./../NURBS/Solid/data_multi_part
      
%  for i=1:size(B,1)
%   for j=1:size(B,2)
%         B{i,j,3}(2)=B{i,j,3}(2)+0.5;
%   end
% end   
%}

pro3;

% load Pso BestSol

deg.p = p; deg.q=q; deg.r = r; clear p q r


[ X_ ] = getSubDivKVValues( Xi, 1); 
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,1);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
 
KV0.Xi=Xi; KV0.Eta=Eta; KV0.Zeta=Zeta;

%%%---------Modification-Optimization design
  A=0;
     a=-0.1;b=8.5;
%    a=-0.25;b=18.5;
for i=(size(B,1)-11):size(B,1)
     A=A+1;
    for j=1:size(B,2)
     if  B{i,j,1}(2)>=10e-10
         B{i,j,1}(2)=a* B{i,j,1}(1)+b;
         B{i,j,2}(2)=B{i,j,1}(2)-0.1;
         B{i,j,1}(4)=1;
     elseif B{i,j,1}(2)<=-10e-10
         B{i,j,1}(2)=-a* B{i,j,1}(1)-b;
         B{i,j,2}(2)=B{i,j,1}(2)+0.1;
         B{i,j,1}(4)=1;
     end
    end       
end
%%% 

% %%%
%   AA=0;aa=[1.1,1.5];
% for i=(size(B,1)-3):2:(size(B,1)-1)
%      AA=AA+1;
%     for j=1:size(B,2)
%      if  abs(B{i,j,1}(2))<=10e-10
%          B{i,j,1}(3)=aa(AA)* B{i,j,1}(3);
%           B{i,j-1,1}(3)=aa(AA)* B{i,j-1,1}(3);
%            B{i,j+1,1}(3)=aa(AA)* B{i,j+1,1}(3);
%      end
%     end       
% end
% %%%

BJ=cell2mat(B);
  n1=size(B,1);m1=size(B,2);l1=size(B,3);
  [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );
B0=B;
B1=B(:,2:size(B,2),:);
% Number of control points after refinement
n = size(B1,1);
m = size(B1,2);
l = size(B1,3);
% B1=B(:,2:m,:);
%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m1-deg.q) * (l-deg.r); % number of elements

nnp = n*m*l; % number of global basis functions

nen = (deg.p+1)*(deg.q+1)*(deg.r+1); % number of local basis functions
ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
% [INN,IEN,AA] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
 [INN,IEN,INB] = BLDINCIEN( deg,n,m,l );
% for i=1:numel(IEN)
%     ie=IEN(i);
%    if (ie>=(n*m1+1)&&ie<=n*m)||(ie>=(n*(2*m-1)+1)&&ie<=n*m*2)||(ie>=(n*(l*m-1)+1)&&ie<=n*m*l)
%           IEN(i)=IEN(i-n*m1);
%    end
%   
% end
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end

%% Material parameters:
E_Y = 210e4;
nu_P = 0.3;
denisty=785;
damp=0.1;
D=hooke_strain(5,E_Y,nu_P);

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
% K = zeros(ndof);
F = zeros(ndof,1);
for e = 1 : nel
    
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk))
        continue
    end
    
    Ke = zeros(nen*3);
    Fe = zeros(nen*3,1);
%     Me = zeros(nen*3);
%     Ce = zeros(nen*3);
%     
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,e,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_,kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                Ke = Ke + Ke_;
%                 
%                 % Build Me
%                  Ni=eye(3);
%                  for s=1:length(R)
%                      shamtx(:,(s-1)*3+1:s*3)=R(s)*Ni;
%                  end
%                  Me_=denisty*shamtx'*shamtx*Jmod;
%                  Ce_=damp*shamtx'*shamtx*Jmod;
%                  Me = Me + Me_;
%                  Ce = Me + Me_;
            end
        end
    end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
%     M(idx,idx) = M(idx,idx) + Me;
%     C(idx,idx) = C(idx,idx) + Ce;
%     F(idx) = F(idx)+Fe;
end

% Apply load at two corner nodes
%     loc_z=0;loc_x=0;   %data_10的
%     loc_y=0;loc_x=0;
%   loc_y=8;loc_x=3;  %chilun
  loc_y=2.5;loc_x=30;  %guan
constNod1 = [];
for i = 9:10
% %        if B{i}(3) == loc_z
% %            if (B{i}(2) == 0.495||B{i}(2) == 0.505)
%           if B1{i}(2) == loc_y
%              if B1{i}(1) == loc_x
%                 constNod=[constNod i];
%              end
%           end
%         
% %         end
  for j=1:size(B1,2)
     
   if abs(B1{i,j,1}(3))<=1e-10
      
    constNod1=[constNod1 INB(i,j,1)];
   end
  end
end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:

F(ID(2,constNod1)) = 3e3;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
%   loc_z=-1.2;
%   loc_y=-1.3;
   loc_x=73;
constNod = [];
for i = 1 : numel(B1)
    if B1{i}(1) == loc_x
        constNod=[constNod i];
    end
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];



%% Solve system
[a,K]=solveq(K,F,bc);

%%
ai=zeros(3*nen,nel);num_node=zeros(nnp,1);
for sq = 1:nel
     ai(:,sq)=a(LM(:,sq)) ;
    
      ni = INN(IEN(1,sq),1); nj = INN(IEN(1,sq),2);nk = INN(IEN(1,sq),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk))
        continue
    end
    gp=1;np=1;
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,sq,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                [ Ke_, kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                strain(:,gp,sq)=kinmtx*ai(:,sq);
                siga=D*strain(:,gp,sq);
                
                stress(gp,sq,1:3)=siga(1:3);
          % von Mises stress
                stress(gp,sq,4)  = sqrt(siga(1)^2+siga(2)^2-...
                                   siga(1)*siga(2)+3*siga(3)^2); 
                stress_n(IEN(gp,sq),:)=stress(gp,sq,:);
                num_node(IEN(gp,sq))= num_node(IEN(gp,sq))+1;
                gp=gp+1;
            end
        end
    end   
end

for i=1:n*m*l
  h=INN(i,:)  ;
  gcoord(i,:)=B1{h(1),h(2),h(3)}(1:3);
%   gcoord(i,2:3)=gcoord(i,2:3).*B1{h(1),h(2),h(3)}(4);
end
for i=1:nnp
    num=num_node(i);
    stress_n(i,:)=stress_n(i,:)./num;
end
% for j=1:l
%  for i=(1+(j-1)*m*n):((j-1)*m*n+n)
%     stress_n(i,:)= (stress_n(i,:)+ stress_n(i+n*(m-1),:))/2;
%     stress_n(i+n*(m-1),:)=stress_n(i,:);
%  end
% end



%%------------环形模型需增加的的应力操作---------

%  mise=zeros(n1*m1*l1,4);ip=size(B1,2)*size(B1,1);
%  for i=1:ip
%      mise((l1*(i-1)+1):(l1*i-1),:)=stress_n(((l1-1)*(i-1)+1):(l1-1)*i,:);
%      mise(l1*i,:)=stress_n((l1-1)*(i-1)+1,:);
%  end
%}

%  mise=zeros(n*m1*l,4);ip=size(B,2)*size(B,1);
%  
%  mise(1:n*m*(l),:)=stress_n;
%  mise((n*m*l)+1):end,:)=stress_n(1:ip,:);
 

%% TRANSFORM the matrix of element(p+1)*(q+1)(r+1) into 2*2*2
    e=0;bi=0;
    p= deg.p ; q=deg.q; r=deg.r ;
for k = 1 : (l-1)
        for j = 1 : (m)
            for i = 1 : (n-1)
                e=e+1;
               for kloc = 0 : 1   %1 ,  3
                   for jloc = 0 : 1   %2 ,  12
                           for iloc = 0 : 1    %3 ,  55                     
                              if (jloc+j)==(m+1)
                               
                                B_j=INB((iloc+i), 1, kloc+k); 
                               else
                                B_j=INB((iloc+i), jloc+j, kloc+k);   % global function number
                              end
                              
                              b = (1+1)*(1+1)*(1+1)-kloc*(1+1)*(1+1) - jloc*(1+1)- iloc ; % local function number
                              BR(b,e) = B_j; % assign connectivity
                           end
                   end
               end
             
            end
        end
end
    BR1=BR;
    BR1(1,:)=BR(2,:);
    BR1(2,:)=BR(1,:);
    BR1(5,:)=BR(6,:);
    BR1(6,:)=BR(5,:);
    
%   tecplot(BR1,gcoord,mise);
     tecplot(BR1,gcoord,stress_n);
%}
%}

%%
u = cell(size(B1)); 
comb = cell(size(B1));
for i = 1 : size(ID,2) %
u{i} = [a(ID(:,i)); 0]; 
comb{i} = B1{i} + u{i}; 
end
       comb(:,m1,:)=comb(:,1,:);
       comb1=cell2mat(comb);


%% plot
plotNurbsSolidElementSimple( KV, B )
title('Initial geometry_the circular plate coupled to a doubly curved barrel shape')   %data_solid10
 plotNurbsSolidElementSimple( KV, comb )
 title('Displacements')
% 
%{
Conn(1:4,:)=IEN(1:4,:);
Conn(5,:)=IEN(8,:);
Conn(6,:)=IEN(7,:);
Conn(7,:)=IEN(6,:);
Conn(8,:)=IEN(5,:);
Conn(9:12,:)=IEN(9:12,:);
Conn(13,:)=IEN(16,:);
Conn(14,:)=IEN(15,:);
Conn(15,:)=IEN(14,:);
Conn(16,:)=IEN(13,:);
Conn(17:20,:)=IEN(17:20,:);
Conn(21,:)=IEN(24,:);
Conn(22,:)=IEN(23,:);
Conn(23,:)=IEN(22,:);
Conn(24,:)=IEN(21,:);
%}
VTKPostProcess(gcoord,BR1','B8','.vtu',stress_n,a);

 disp0=a;K0=K;L=chol(K); ndof0=ndof;
 
  save initial  disp0 ndof0  K0 L B0 deg KV0 D
toc

