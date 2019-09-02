function  [Jss,RR,G]= Jx_s_smatrix(c,zeta,dS_ds,dS_dt,R,v3,dR_dxi,v1,v2,Jss)
  ti=1;  
      c1=norm(c).*(cross(dS_dss, dS_dt)+cross(dS_ds,dS_dst));
      c2=c.*norm(cross(dS_dss, dS_dt)+cross(dS_ds,dS_dst));
      v3_s=(c1-c2)/norm(c)^2;
      c11=norm(c).*(cross(dS_dtt, dS_ds)+cross(dS_dt,dS_dst));
      c22=c.*norm(cross(dS_dtt, dS_ds)+cross(dS_dt,dS_dst));  
      v3_t=(c11-c22)/norm(c)^2;
    e1=[1;0;0];e2=[0;1;0];
       if cross(e2,v3(lo,:))==zeros(3,1)
          v2_s=cross(e1,v3_s);
         v1_s=cross( v2_s,v3)+cross(v3_s,v2);
        else
          v1_s=cross(e2,v3_s);
          v2_s=cross(v1,v3_s)+cross(v1_s,v3);
       end
      
     
      for i=1:length(Jss)
           Jss(3,i)=v3(i)*R(lo)*ti/2;
          
           Jss(1,i)= dS_ds(i)+zeta*(v3_s*R(lo)+v3(i)*dR_dxi(lo,1))*ti/2;
           Jss(2,i)=  dS_dt(i)+zeta*(v3_t*R(lo)+v3(i)*dR_dxi(lo,2))*ti/2;
          
      end
      
       Jsinv = inv(Jss);
       O=zeros(3);
       G=[Jsinv,    O,    O;
           O,    Jsinv ,  O;
           O,       O,   Jsinv];
       RR=zeros(9,5);dd=zeros(3);
       for i=1:length(RR)
           a=round(i/3)+1;b=mod(i/3);
           if b~=0
           RR(i,a)=dR_dxi(lo,1);
           RR(i,4)=-zeta*(v2_s(1)*R(lo)+v2(1)*dR_dxi(lo,1))*ti/2;
           RR(i,5)=zeta*(v1_s(1)*R(lo)+v1(1)*dR_dxi(lo,1))*ti/2;
           else
              RR(i,4)=  v2(i)*R(lo)*ti/2; 
              RR(i,5)=  v1(i)*R(lo)*ti/2;
           end  
       end
          
     