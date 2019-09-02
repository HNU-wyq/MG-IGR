function [objtive,aT]=obj_strewdriver(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1)
%% Load data
addpath(genpath('C:\Users\Administrator\Desktop\��ҵ����ljj\����\calfem-matlab-iga-master-screwdriver'));
count=1;
%{
%  data_l;
%      run ./../NURBS/Solid/data_solid_chilun
%          run ./../NURBS/solid/data_solid10
          run ./../NURBS/Solid/data_part2
%             run ./../NURBS/Solid/data_multi_part
     

pro3;

% load Pso BestSol

deg.p = p; deg.q=q; deg.r = r; clear p q r



[ E_ ] = getSubDivKVValues( Eta,1);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
% [ X_ ] = getSubDivKVValues( Xi, 1); 
% [B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
%}

load initial 
B=B0;
B01=B1;
p=deg.p ; q=deg.q; r=deg.r ; 

%% ---------initial Modification-initial design
 
     a=-0.1;b=9;
%    a=-0.25;b=18.5;
% for i=(size(B,1)-6):size(B,1)
%      Aj=Aj+1;
%     for j=1:size(B,2)
%      if  B{i,j,1}(2)>=10e-10
%          B{i,j,1}(2)=a* B{i,j,1}(1)+b;
%          B{i,j,2}(2)=B{i,j,1}(2)-0.1;
%          B{i,j,1}(4)=1;
%      elseif B{i,j,1}(2)<=-10e-10
%          B{i,j,1}(2)=-a* B{i,j,1}(1)-b;
%          B{i,j,2}(2)=B{i,j,1}(2)+0.1;
%          B{i,j,1}(4)=1;
%      end
%     end       
% end
% for i=(size(B,1)-4):size(B,1)
%      Aj=Aj+1;
%     for j=1:size(B,2)
%          B{i,j,1}(2)=B{size(B,1),j,1}(2);
%     end       
% end

 Aj=0;
 for i=(size(B,1)-6):size(B,1) 
     Aj=Aj+1;
%     for j=1:size(B,2)
%      if  B{i,j,1}(2)>=10e-10
%          B{i,j,1}(2)=2;
%      elseif B{i,j,1}(2)<=-10e-10
%          B{i,j,1}(2)=-2;
%      end
%     end 
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
    mid=cell2mat(B(i,:,1));
    nub=find(abs(abs(mid(3,:))-max(mid(3,:)))<1e-8);
    for k=1:numel(nub)
     if  B{i,nub(k),1}(2)>=10e-10    
          B{i,nub(k),1}(2)=B{i,nub(k),1}(2)-x(1);
          B{i,nub(k),2}(2)=B{i,nub(k),2}(2)-x(1);             %��ͷ���Ż�Ϊ͹���Σ�8���Σ� һ����Ʊ����ı�6����Ƶ� ��0.8-1.8��
     elseif B{i,nub(k),1}(2)<=-10e-10
          B{i,nub(k),1}(2)=B{i,nub(k),1}(2)+x(1);
          B{i,nub(k),2}(2)=B{i,nub(k),2}(2)+x(1);
     end
    end
end

%% optimization the handle of screwdriver
% x=[0.9,1.10,1.15,0.95,1.15,1.1,1.25,0.85,1.35,1.1,0.7,0.75];
Ajs=2;
for i=(size(B,1)-12):(size(B,1)-7)  %���ˣ�������ͷ�����Ż�  6����Ƶ㣬6����Ʊ��� ��0.8-1.2��
    
    for j=1:size(B,2)
        
         B{i,j,1}(2:3)=x(Ajs).* B{i,j,1}(2:3);
    end    
      Ajs=Ajs+1;
end
for i=12:14   %�����ڰ���   3����Ʊ�����3����Ƶ�  ��0.5,1.5��
     
    for j=1:size(B,2)
        
         B{i,j,1}(2:3)=x(Ajs).* B{i,j,1}(2:3);    
    end    
     Ajs=Ajs+1;
end
% Ajs
% x=[1.26,2.29];
%  ni=(size(B,1)-2);
%     for j=1:size(B,2)
%      if  B{ni,j,1}(2)<=10e-10 && B{ni,j,1}(2)>0
%          B{ni,j,1}(3)=x(1)* B{ni,j,1}(3);
%           B{ni,j-1,1}(3)=x(1)* B{ni,j-1,1}(3);
%            B{ni,j+1,1}(3)=x(1)* B{ni,j+1,1}(3);
%      elseif B{ni,j,1}(2)>=-10e-10 && B{ni,j,1}(2)<0
%          B{ni,j,1}(3)=x(1)* B{ni,j,1}(3);
%           B{ni,j-1,1}(3)=x(1)* B{ni,j-1,1}(3);
%            B{ni,j+1,1}(3)=x(1)* B{ni,j+1,1}(3);
%          
%      end
%     end 
%      ni=ni+1;
%     for j=1:size(B,2)
%      if  B{ni,j,1}(2)<=10e-10 && B{ni,j,1}(2)>0
%          B{ni,j,1}(3)=x(2)* B{ni,j,1}(3);
%           B{ni,j-1,1}(3)=x(2)* B{ni,j-1,1}(3);
%            B{ni,j+1,1}(3)=x(2)* B{ni,j+1,1}(3);
%      elseif B{ni,j,1}(2)>=-10e-10 && B{ni,j,1}(2)<0
%          B{ni,j,1}(3)=x(2)* B{ni,j,1}(3);
%           B{ni,j-1,1}(3)=x(2)* B{ni,j-1,1}(3);
%            B{ni,j+1,1}(3)=x(2)* B{ni,j+1,1}(3);        
%      end
%     end 
    
B0=B;Xi=KV0.Xi; Eta=KV0.Eta; Zeta=KV0.Zeta;
[ X_ ] = getSubDivKVValues( Xi, 1); 
[B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
%%%

BJ=cell2mat(B);
  n1=size(B,1);m1=size(B,2);l1=size(B,3);
  [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );

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
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) ��ϵͳ�е�Ԫ�ڵ����˵�����ɶȱ�š�
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
% K = zeros(ndof); % Needs to be changed to sparse for large problems!!
K = sparse(ndof,ndof);
F = zeros(ndof,1);
for e = 1 : nel
     je=1;
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk)))&&((e<((n-p)*(m-q)+1)||e>((n-p)*(m-q+1)))&&(e<((n-p)*(m-q)+1+(n-p)*(m-q+1))||e>((n-p)*(m-q+1)*2)))
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
                [ R,dR_dx,Jdet ] = Shape_function_repoint( GP,e,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                J(je,e)=Jmod;
                % Build Ke
                [ Ke_,kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                Ke = Ke + Ke_;
                Bmatrix{je,e}=kinmtx;
                je=je+1; 
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
%     tic
ke=sparse(repmat(idx,1,numel(idx)),kron(idx,ones(1,numel(idx))),reshape(Ke,1,numel(Ke)),ndof,ndof);
%     K(idx,idx) = K(idx,idx) + Ke;
    K = K + ke;
%     toc
%     M(idx,idx) = M(idx,idx) + Me;
%     C(idx,idx) = C(idx,idx) + Ce;
%     F(idx) = F(idx)+Fe;
end

% Apply load at two corner nodes
%     loc_z=0;loc_x=0;   %data_10��
%     loc_y=0;loc_x=0;
%   loc_y=8;loc_x=3;  %chilun
  loc_y=2.5;loc_x=30;  %guan
constNod1 = [];
for i = 6:1:16
% %        if B{i}(3) == loc_z
% %            if (B{i}(2) == 0.495||B{i}(2) == 0.505)
%           if B1{i}(2) == loc_y
%              if B1{i}(1) == loc_x
%                 constNod=[constNod i];
%              end
%           end
%         
% %         end
%    if abs(B1{i,j,1}(2))<=1e-10  
  for j=1:size(B1,2)
    
    constNod1=[constNod1 INB(i,j,1)];
    
    ay=B1{INB(i,j,1)}(2);az=B1{INB(i,j,1)}(3);
    F(ID(2,INB(i,j,1)) )=-az/sqrt(ay^2+az^2)*2000;     %ת�ص�Ч�ڵ�������ʱ����ת��P
    F(ID(3,INB(i,j,1)) )=ay/sqrt(ay^2+az^2)*2000;      %ת�ص�Ч�ڵ�������ʱ����ת��P
  end

end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:
% F(ID(2,constNod1) )=- 3e2;
% F(ID(2,constNod1(1:2:end))) =- 3e2;
% F(ID(2,constNod1(2:2:end))) = 3e2;


%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
%   loc_z=-1.2;
%   loc_y=-1.3;
 loc_x=73;
 constNod = [];
% for i = 1 : numel(B1)
%     if B1{i}(1) == loc_x
%         constNod=[constNod i];
%     end
% end

for i=(size(B,1)-6):size(B,1)
    for j=1:size(B1,2)
%         for k=1:size(B1,3)
        if B1{i,j,1}(2)>0.01&&B1{i,j,1}(3)<-0.01
           constNod=[constNod INB(i,j,1)]; 
        elseif  B1{i,j,1}(2)<-0.01&&B1{i,j,1}(3)>0.01
           constNod=[constNod INB(i,j,1)]; 
        end
    end
%      constNod=constNod(1:2:end);    %��͹�����ʱ
        %end
%     mid=cell2mat(B1(i,:,1));
%     nub=find(abs(abs(mid(3,:))-max(mid(3,:)))<1e-8);
%     for k=1:numel(nub)
%      if  B{i,nub(k),1}(2)>=10e-10
%         constNod=[constNod INB(i,nub(k),1)]; 
%      elseif B{i,nub(k),1}(2)<=10e-10
%         constNod=[constNod INB(i,nub(k),1)]; 
%      end
%     end
end

if x(1)==0
      constNod=constNod(1:2:end);    %��͹�����ʱ  
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];


%% Solve system
% [a,K]=solveq(K,F,bc);

nn=length(bc1);
for i3=1:nn
    c=bc1(i3);
    K(c,:)=0;
    K(:,c)=0;
    K(c,c)=1; 
    F(bc1(i3))=0;
end
% K=sparse(K);
%% Solve system
% [a,K]=solveq(K,F,bc);

[aT,Kk1, Pp2] = Test22(K,F,P2,disp0,B1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
clear Kk1 Pp1 Pp2   P2 INBz1 INNz1 Br1 Bz11 Bz  nsub ssub idsub


%%
ai=zeros(3*nen,nel);num_node=zeros(nnp,1);
for sq = 1:nel
     ai(:,sq)=aT(LM(:,sq)) ;
    
      ni = INN(IEN(1,sq),1); nj = INN(IEN(1,sq),2);nk = INN(IEN(1,sq),3);
    
    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk)))&&((sq<((n-p)*(m-q)+1)||sq>((n-p)*(m-q+1)))&&(sq<((n-p)*(m-q)+1+(n-p)*(m-q+1))||sq>((n-p)*(m-q+1)*2)))
        continue
    end
  
   for gp=1:size(IEN,1)
         %{
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
%                 [ R,dR_dx,Jdet ] = Shape_function( GP,sq,deg,B1,KV,INN,IEN);
                 [ R,dR_dx,Jdet ] = Shape_function_repoint( GP,sq,deg,B1,KV,INN,IEN);
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                [ Ke_, kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                %}
                kinmtx=Bmatrix{gp,sq};
                strain(:,gp,sq)=kinmtx*ai(:,sq);
                siga=D*strain(:,gp,sq);
        
                stress(gp,sq,1:3)=siga(1:3);
          % von Mises stress
                stress(gp,sq,4)  = sqrt(((siga(1)-siga(2))^2+(siga(3)-siga(2))^2+...
                                   (siga(1)-siga(3))^2)/2); 
                stress_n(IEN(gp,sq),:)=stress(gp,sq,:);
%                 num_node(IEN(gp,sq))= num_node(IEN(gp,sq))+1;
%                 gp=gp+1;
%             end
%         end
%     end   
    end
end

% for i=1:n*m*l
%   h=INN(i,:)  ;
% gcoord(i,:)=B1{h(1),h(2),h(3)}(1:3);
% % gcoord(i,2:3)=B1{h(1),h(2),h(3)}(2:3).*B1{h(1),h(2),h(3)}(4);
% end

for jj=1:size(IEN,1)
   for ii=1:nel
        num_node(IEN(jj,ii))= num_node(IEN(jj,ii))+1;
    end
end
 
% for i=1:nnp
%     num=num_node(i);
%     stress_n(i,:)=stress_n(i,:)./num;
% end
stress_n=bsxfun( @rdivide,stress_n,num_node);

objtive=max(stress_n(:,4));
% for j=1:l
%  for i=(1+(j-1)*m*n):((j-1)*m*n+n)
%     stress_n(i,:)= (stress_n(i,:)+ stress_n(i+n*(m-1),:))/2;
%     stress_n(i+n*(m-1),:)=stress_n(i,:);
%  end
% end



%%------------����ģ�������ӵĵ�Ӧ������---------

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

%{


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
     tecplot(BR1,gcoord,stress_n,a);
%}
%}

% %
% u = cell(size(B1)); 
% comb = cell(size(B1));
% for i = 1 : size(ID,2) %
% u{i} = [a(ID(:,i)); 0]; 
% comb{i} = B1{i} + u{i}; 
% end
%        comb(:,m1,:)=comb(:,1,:);
%        comb1=cell2mat(comb);
% 
% 
% %% plot
% plotNurbsSolidElementSimple( KV, B )
% title('Initial geometry_the circular plate coupled to a doubly curved barrel shape')   %data_solid10
%  plotNurbsSolidElementSimple( KV, comb )
%  title('Displacements')
% % 
%  disp0=a;K0=K;L=chol(K); ndof0=ndof;
%  
%   save initial  disp0 ndof0  K0 L B0 deg KV0 D
% toc
%}

