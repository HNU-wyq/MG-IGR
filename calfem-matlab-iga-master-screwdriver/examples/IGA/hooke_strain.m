function [D]=hooke_strain(iopt,E_Y,nu_P)

%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%
%  Synopsis:
%     [matmtrx]=fematiso(iopt,elastic,poisson) 
%
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - ���ο�Ԫ
%     iopt=5 - three dimensional analysis
%------------------------------------------------------------------------

 if iopt==1        % plane stress
   D= E_Y/(1-nu_P*nu_P)* ... 
   [1  nu_P 0; ...
   nu_P  1  0; ...
   0  0  (1-nu_P)/2];

 elseif iopt==2     % plane strain
   D= E_Y/((1+nu_P)*(1-2*nu_P))* ...
   [(1-nu_P)  nu_P 0; 
   nu_P  (1-nu_P)  0;
   0  0  (1-2*nu_P)/2];

 elseif iopt==3     % axisymmetry��Գ�
   D= E_Y/((1+nu_P)*(1-2*nu_P))* ...
   [(1-nu_P)  nu_P  nu_P  0; 
   nu_P  (1-nu_P)   nu_P  0;
   nu_P  nu_P  (1-nu_P)   0;
   0    0    0   (1-2*nu_P)/2];
  elseif iopt==4   % ���ο�Ԫ
     k=1.2;
     D= E_Y/(1-nu_P*nu_P)* ...
   [1  nu_P 0  0  0  0; ...
   nu_P  1  0  0  0  0; ...
   0     0  0  0  0  0;
   0  0  0  (1-nu_P)/2 0  0;...
   0  0  0  0 (1-nu_P)/2*k  0;...
   0  0  0  0  0  (1-nu_P)/2*k];
 else     % three-dimension
   D= E_Y/((1+nu_P)*(1-2*nu_P))* ...
   [(1-nu_P)  nu_P  nu_P   0   0    0; 
   nu_P  (1-nu_P)   nu_P   0   0    0;
   nu_P  nu_P  (1-nu_P)    0   0    0;
   0    0    0    (1-2*nu_P)/2   0    0;
   0    0    0    0    (1-2*nu_P)/2   0;
   0    0    0    0    0   (1-2*nu_P)/2];

 end
