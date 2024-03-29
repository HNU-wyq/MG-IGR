function [ Ke,B_m,dme] = shell_Build_K_Local( dR_dx,Jmod,N1,ndof,Zt,thetai,Qmb,nl,Gtw,Glw,h,density)% Generate B matrix
% kinmtx=zeros(6,nen*3);
% Ke=zeros(5*ndof);
% %% membrane + dending 
% Qmb = [El/(1-vlt*vtl) vtl*El/(1-vlt*vtl) 0;
%     vtl*El/(1-vlt*vtl) Et/(1-vlt*vtl)  0;
%     0     0         Glt];
%[B] matrix bending

B_b=zeros(3,5*ndof);
B_m=zeros(3,5*ndof);
% B_b(1,ndof+1:2*ndof) = dR_dx(:,1)';
% B_b(2,2*ndof+1:3*ndof) =dR_dx(:,2)';
% B_b(3,ndof+1:2*ndof) = dR_dx(:,2)';
% B_b(3,2*ndof+1:3*ndof) = dR_dx(:,1)';

B_b(1,4:5:5*ndof) = dR_dx(:,1)';
B_b(2,5:5:5*ndof) =dR_dx(:,2)';
B_b(3,4:5:5*ndof) = dR_dx(:,2)';
B_b(3,5:5:5*ndof) = dR_dx(:,1)';
% [B] matrix membrane

% B_m(1,3*ndof+1:4*ndof) = dR_dx(:,1)';  %ԭ����
% B_m(2,4*ndof+1:5*ndof) = dR_dx(:,2)';
% B_m(3,3*ndof+1:4*ndof) = dR_dx(:,2)';
% B_m(3,4*ndof+1:5*ndof) = dR_dx(:,1)';

B_m(1,1:5:5*ndof) = dR_dx(:,1)';
B_m(2,2:5:5*ndof) = dR_dx(:,2)';
B_m(3,1:5:5*ndof) = dR_dx(:,2)';
B_m(3,2:5:5*ndof) = dR_dx(:,1)';


th = thetai;
kmm=zeros(5*ndof);
kmb=zeros(5*ndof);
kbm=zeros(5*ndof);
kff=zeros(5*ndof);
kcc=zeros(5*ndof);   %shear
for i = 1:nl    %layers
    thi = th(i);   %angles
%     T = [cos(thi)^2 sin(thi)^2 2*cos(thi)*sin(thi);
%         sin(thi)^2 cos(thi)^2 -2*cos(thi)*sin(thi);
%         -cos(thi)*sin(thi) cos(thi)*sin(thi) cos(thi)^2-sin(thi)^2];
    Ti=[cos(thi)^2 sin(thi)^2 -2*cos(thi)*sin(thi);
        sin(thi)^2 cos(thi)^2 2*cos(thi)*sin(thi);
        cos(thi)*sin(thi) -cos(thi)*sin(thi) cos(thi)^2-sin(thi)^2];
    
    z1 = Zt(i,1);   %Zt
    z2 = Zt(i,2);
    Qb=Ti*Qmb*Ti.';

    kmm =kmm+B_m'*Qb*B_m*det(Jmod)*(z2-z1);
    kmb =kmb+B_m'*Qb*B_b*det(Jmod)*(z2^2-z1^2)/2;
    kbm =kbm+B_b'*Qb*B_m*det(Jmod)*(z2^2-z1^2)/2;
    kff =kff+B_b'*Qb*B_b*det(Jmod)*(z2^3-z1^3)/3;

  
end
kmbf=kmm+kmb+kbm+kff;

%% shear 
Bc=zeros(2,5*ndof);
% Bc(1,1:ndof)=dR_dx(1,:)';
% Bc(2,1:ndof)=dR_dx(2,:)';
% Bc(1,ndof+1:2*ndof)=N1;
% Bc(2,2*ndof+1:3*ndof)=N1;
Bc(1,3:5:5*ndof)=dR_dx(:,1)';
Bc(2,3:5:5*ndof)=dR_dx(:,2)';
Bc(1,4:5:5*ndof)=N1;
Bc(2,5:5:5*ndof)=N1;

kapa=pi*pi/12;
% kapa=5/6;
cs=kapa*[Gtw 0;
         0 Glw];
th = thetai;
kcc=zeros(5*ndof);

for i = 1:nl
    thi = th(i);

    Tt=[cos(thi) sin(thi);
        -sin(thi) cos(thi)];
    Dc=Tt*cs*Tt';

    kcc =kcc+Bc'*Dc*Bc*det(Jmod)*(Zt(i,2)-Zt(i,1));
end

Ke = kmbf+kcc;
% %stiffness matrix
% Ke = kinmtx'*D*kinmtx*Jmod;

I=h^3/12;
me1 = N1*N1'*det(Jmod)*density*h;
me2 = N1*N1'*det(Jmod)*density*I;
dme=zeros(5*ndof);
a=1:5:5*ndof;
dme(a,a)=dme(a,a)+me1;
dme(a+1,a+1)=dme(a+1,a+1)+me1;
dme(a+2,a+2)=dme(a+2,a+2)+me1;
dme(a+3,a+3)=dme(a+3,a+3)+me2;
dme(a+4,a+4)=dme(a+4,a+4)+me2;

end