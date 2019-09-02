function [ Ke ] = Build_K_Local_shell( dR_dx,Jmod,D,nen )
% PURPOSE:
% Calculate the contribution from the current integration point 
% to the linear elastic stiffness matrix for a solid 
% NURBS element.
%
% INPUT: dR_dx = vector of basis function derivatives (nen x 3)
%
%        Jmod  = Jacobian determinant
%
%        D     = constitutive matrix (6 x 6)
%
%        nen   = number of local basis functions
%
% OUTPUT: Ke = element stiffness matrix contribution (nen*3 x nen*3)

%————壳单元分两部分Bv与Bw
% Generate B matrix
Bv=zeros(6,nen*3);
Bv(1,1:3:end) = dR_dx(:,1);
Bv(2,2:3:end) = dR_dx(:,2);
Bv(3,3:3:end) = dR_dx(:,3);
Bv(4,1:3:end) = dR_dx(:,2);
Bv(4,2:3:end) = dR_dx(:,1);
Bv(5,2:3:end) = dR_dx(:,3);
Bv(5,3:3:end) = dR_dx(:,2);
Bv(6,1:3:end) = dR_dx(:,3);
Bv(6,3:3:end) = dR_dx(:,1);

Bv
Bw=zeros(6,nen*3);
y_=[0,0,0.1];h=0.1
dR_dx_=dR_dx.*Zeta;
Bw(1,1:3:end) =-dR_dx_(:,2).*y_(3);
Bw(1,3:3:end) =dR_dx_(:,2).*y_(1);
Bw(2,1:3:end) =dR_dx_(:,3).*y_(2);
Bw(2,2:3:end) =-dR_dx_(:,3).*y_(1);
Bw(3,1:3:end) =-dR_dx_(:,1).*y_(3);
Bw(3,2:3:end) =dR_dx_(:,2).*y_(3);
Bw(3,3:3:end) =dR_dx_(:,1).*y_(1)-dR_dx_(:,2).*y_(2);
Bw(4,1:3:end) =dR_dx_(:,2).*y_(2)-dR_dx_(:,3).*y_(3);
Bw(4,2:3:end) =-dR_dx_(:,2).*y_(1);
Bw(4,3:3:end) =dR_dx_(:,3).*y_(1);
Bw(5,1:3:end) =dR_dx_(:,1).*y_(2);
Bw(5,2:3:end) =dR_dx_(:,3).*y_(3)-dR_dx_(:,1).*y_(1);
Bw(5,3:3:end) =-dR_dx_(:,3).*y_(2);

Bw=Bw./2


% Contrubution to element stiffness matrix
B=[Bv,Bw];
Ke = B'*D*B*Jmod;

end