function [ Bm,J,R ] = Shape_function_membrane( GP,e,deg,B,KV,INC,IEN)
% [ R,dR_dx,J ] = Shape_function( GP,e,deg,B,KV,INN,IEN)
%-------------------------------------------------------------
% PURPOSE:
% For a solid NURBS element, calculate the the vector of local 
% shape functions R   and an array of their derivatives dR_dx    at
% the current gauss point, and     the Jacobian determinant J.
% 
% INPUT: The quadrature point GP, containing:
%        GP.Xi_tilde, GP.Eta_tilde, and GP.Zeta_tilde 
%
%        element number e, 
%
%        polynomial orders deg, containing:
%        deg.p, deg.q, and deg.r
%
%        control net B{n,m,l}(x,y,z,w), with a cell structure 
%        with n points in Xi direction, m points in Eta 
%        direction, and l points in Zeta direction. Weights 
%        are stored as the fourth component.
%
%        knot vectors KV, containing the clamped non-uniform 
%        knot vectors: 
%        KV.Xi, KV.Eta, and KV.Zeta
%
%        and connectivity arrays INC and IEN.
% 
% OUTPUT: R = vector of basis functions (nen x 1)
%         dR_dx = vector of basis function derivatives (nen x 3)
%         J = Jacaboian determinant (1)
%-------------------------------------------------------------

p = deg.p; q = deg.q;   r=deg.r; 

% number of local basis functions:
% nen = (p+1)*(q+1)*(r+1); 
  nen = (p+1)*(q+1);                       
% NURBS coordinates; convention consistent with Algorithm 7
ni = INC(IEN(1,e),1);
nj = INC(IEN(1,e),2);
% nk = INC(IEN(1,e),3);

% Calculate parametric coordinates from parent element coordinates
% Knot vectors KV.Xi, KV.Eta, and KV.Zeta and
% parent element coordinates xi_tilde, eta_tilde, zeta_tilde
% are given as input
xi = ((KV.Xi(ni+1)-KV.Xi(ni))*GP.xi_tilde ...
+ (KV.Xi(ni+1)+KV.Xi(ni))) / 2;
eta = ((KV.Eta(nj+1)-KV.Eta(nj))*GP.eta_tilde ...
+ (KV.Eta(nj+1)+KV.Eta(nj))) / 2;
zeta = GP.zeta_tilde;   %超参壳元的gamma
% zeta = ((KV.Zeta(nk+1)-KV.Zeta(nk))*GP.zeta_tilde ...
% + (KV.Zeta(nk+1)+KV.Zeta(nk))) / 2;

% Calculate univariate B-spline functions using (2.1) and (2.2)
% and their derivatives using (2.12)
N1 = Der1BasisFun(ni-1, xi,p,KV.Xi)'; % xi-dir.
N2 = Der1BasisFun(nj-1,eta,q,KV.Eta)'; % eta-dir.
     %N3 = Der1BasisFun(nk-1,zeta,r,KV.Zeta)'; % zeta-dir.
N =  N1(:,1);    %三维的形函数N,M,L
dN_dxi = N1(:,2);
M = N2(:,1);
dM_deta = N2(:,2);
   %L = N3(:,1);
   %dL_dzeta = N3(:,2);
clear N1 N2 

% Build numerators and denominators (in local numbering)
x = zeros(1,nen); y = zeros(1,nen);z = zeros(1,nen);
R = zeros(nen,1);   % Array of trivariate NURBS basis functions
dR_dxi = zeros(nen,2); % Trivariate NURBS function derivatives
                       % w.r.t. parametric coordinates
dR_dxdeta=   zeros(nen,1);                    

loc_num = 0; % Local basis function counter


%-------求Ni,s即形函（the rational shape function）的参数s,t的一阶导
  
    for j = 0 : q
        for i = 0 : p
            loc_num = loc_num + 1;
            
%            R(loc_num) = N(p+1-i)*M(q+1-j)*L(r+1-k) ...
%                         * B{ni-i,nj-j,nk-k}(4); % Function numerator (N*M*L*w)
            R(loc_num) = N(p+1-i)*M(q+1-j) ...
                         * B{ni-i,nj-j}(4);
            % Get coordinates in local numbering
            x(loc_num) = B{ni-i,nj-j}(1) ;
            y(loc_num) = B{ni-i,nj-j}(2);
            z(loc_num) = B{ni-i,nj-j}(3);
 
             dR_dxi(loc_num,1) = dN_dxi(p+1-i)*M(q+1-j) ...
            * B{ni-i,nj-j}(4); % Derivative numerator (dN*M*L*w)  B{ni-i,nj-j}(4)这部分是权函数
        
             dR_dxdeta(loc_num) = dN_dxi(p+1-i)*dM_deta(q+1-j) ...
            * B{ni-i,nj-j}(4); 
             
             dR_dxi(loc_num,2) = N(p+1-i)*dM_deta(q+1-j)...
            * B{ni-i,nj-j}(4); % Derivative numerator (N*dM*L*w)

     

            %    dR_dxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dL_dzeta(r+1-k) ...
               %   * B{ni-i,nj-j,nk-k}(4); % Derivative numerator (N*M*dL*w)
        end
    end

W = sum(R); % Function denominator (Sum(N*M*L*w))
dW_dxi = sum(dR_dxi(:,1)) % Derivative denominator (Sum(dN*M*L*w))
dW_deta = sum(dR_dxi(:,2)) % Derivative denominator (Sum(N*dM*L*w))
%dW_dzeta = sum(dR_dxi(:,3)); % Derivative denominator (Sum(N*M*dL*w))
dW_dxdeta = sum(dR_dxdeta(:,1));  
            
% Divide by denominators to complete definitions of functions
% and derivatives w.r.t. parametric coordinates
dR_dxi(:,1) = (dR_dxi(:,1)*W - R*dW_dxi) / W^2   ;      %s
dR_dxi(:,2) = (dR_dxi(:,2)*W - R*dW_deta) / W^2 ;        %t
%dR_dxi(:,3) = (dR_dxi(:,3)*W - R*dW_dzeta) / W^2;

R = R/W;


%%-----------------------5个自由度的mindlin壳---------------------------------------------------------
Bm=zeros(6,nen*5);ti=0.1;

v3=zeros(3,1);v1=zeros(3,1);v2=zeros(3,1);
dS_ds=zeros(3,1); dS_dt=zeros(3,1);
%----求方向余弦及其对应一阶导
e1=[1;0;0]; e2=[0;1;0]; Jss=zeros(3);

      dS_ds=[x;y;z] * dR_dxi(:,1);
      dS_dt=[x;y;z] * dR_dxi(:,2);
      c= cross(dS_ds,dS_dt);
      v3=c/norm(c);
        if cross(e2,v3)==zeros(3,1)
             v2=cross(v3,e1);
             v1=cross(v3,v2);
        else
             v1=cross(e2,v3);
             v2=cross(v3,v1);
        end
        
[H,T]=transmatrix(v1,v2,v3);



 %--------
dS_dst=zeros(3,1);dS_dtt=zeros(3,1);dS_dss=zeros(3,1); 
dS_ss=-2*W*dW_dxi*(dR_dxi(:,1)*W - R*dW_dxi)/ W^4; 
front=(dR_dxdeta*W+dW_deta*dR_dxi(:,1)-dR_dxi(:,2)*dW_dxi-R* dW_dxdeta)*W^2;
dS_st=(front-2*W*dW_deta*(dR_dxi(:,1)*W - R*dW_dxi))  / W^4;         
dS_tt=-2*W*dW_deta*(dR_dxi(:,1)*W - R*dW_dxi)/ W^4; 
  
dS_dss=[x;y;z]*dS_ss;
dS_dst=[x;y;z]*dS_st;
dS_dtt=[x;y;z]*dS_tt;

[v3_s,v3_t,v2_s,v1_s,v2_t,v1_t]=stma(dS_dtt,dS_dst,dS_dss,c,v1,v2,v3,dS_dt,dS_ds);

      dR_ds=sum(dR_dxi(:,1));  dR_dt=sum(dR_dxi(:,2));
      for i=1:length(Jss)

           Jss(1,i)= dS_ds(i)+zeta*(v3_s(i)*sum(R)+v3(i)*dR_ds)*ti/2;
           Jss(2,i)=  dS_dt(i)+zeta*(v3_t(i)*sum(R)+v3(i)*dR_dt)*ti/2;
           Jss(3,i)=v3(i)*sum(R)*ti/2;
      end
      
       Jsinv = inv(Jss);
       O=zeros(3);
       G=[Jsinv, O, O; O, Jsinv ,O; O, O,Jsinv];
       
for lo=1:loc_num
   
      RR=zeros(9,5);
       for i=1:length(RR)
           a=floor(i/3)+1;
           b=mod(i,3);
           
               if b==1
           RR(i,a)=dR_dxi(lo,1);
           RR(i,5)=-zeta*(v2_s(a)*sum(R)+v2(a)* dR_ds)*ti/2;
           RR(i,4)= zeta*(v1_s(a)*sum(R)+v1(a)* dR_ds)*ti/2;
               elseif b==2
           RR(i,a)=dR_dxi(lo,2);
           RR(i,5)=-zeta*(v2_t(a)*sum(R)+v2(a)* dR_dt)*ti/2;
           RR(i,4)= zeta*(v1_t(a)*sum(R)+v1(a)*dR_dt)*ti/2; 
               
               elseif b==0
              RR(i,4)= -v2(a-1)*sum(R)*ti/2; 
              RR(i,5)=  v1(a-1)*sum(R)*ti/2;
               end  
       end

       BB=T*H*G*RR;
       Bm(1:6,5*(lo-1)+1:5*lo)=BB;

            
end
 JR=[0.5*(KV.Xi(ni+1)-KV.Xi(ni)),  0, 0;
     0,  0.5*(KV.Eta(nj+1)-KV.Eta(nj)),0;
     0,  0,                           1];
 J=JR*Jss;
 
%{       
% Gradient of mapping from parameter space to physical space
dx_dxi = [x;y;z] * dR_dxi;

% Compute derivatives of basis functions
% the derivate of R with respect to physical coordinates x,y,z
dR_dx = dR_dxi/dx_dxi;     %---dR_dxi * inv(dx_dxi)

% Gradient of mapping from parent element to parameter space
dxi_dtildexi=zeros(3); % Derivative of parametric coordinates
                       % w.r.t. parent element coordinates
dxi_dtildexi(1,1) = (KV.Xi(ni+1)-KV.Xi(ni))/2;
dxi_dtildexi(2,2) = (KV.Eta(nj+1)-KV.Eta(nj))/2;
dxi_dtildexi(3,3) = (KV.Zeta(nk+1)-KV.Zeta(nk))/2;

% Compute the Jacobian
J_mat = dx_dxi*dxi_dtildexi;  %(3*3)

% Compute Jacobian determinant
J = det(J_mat);
       %}

end