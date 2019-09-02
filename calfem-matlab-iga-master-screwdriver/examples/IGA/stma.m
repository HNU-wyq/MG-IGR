function [v3_s,v3_t,v2_s,v1_s,v2_t,v1_t]=stma(dS_dtt,dS_dst,dS_dss,c,v1,v2,v3,dS_dt,dS_ds)



   ti=1;  
      c1=norm(c).*(cross(dS_dss, dS_dt)+cross(dS_ds,dS_dst));
      c2=c.*norm(cross(dS_dss, dS_dt)+cross(dS_ds,dS_dst));
      v3_s=(c1-c2)/norm(c)^2;
      c11=norm(c).*(cross(dS_dtt, dS_ds)+cross(dS_dt,dS_dst));
      c22=c.*norm(cross(dS_dtt, dS_ds)+cross(dS_dt,dS_dst));  
      v3_t=(c11-c22)/norm(c)^2;
      e1=[1;0;0];e2=[0;1;0];
       if cross(e2,v3_s)==zeros(3,1)
          v2_s=cross(v3_s,e1);
          v1_s=cross( v2_s,v3)+cross(v2,v3_s);
        else
          v1_s=cross(e2,v3_s);
          v2_s=cross(v3_s,v1)+cross(v3,v1_s);
       end
       
       if cross(e2,v3_t)==zeros(3,1)
          v2_t=cross(v3_t,e1);
          v1_t=cross( v2_t,v3)+cross(v2,v3_t);
        else
          v1_t=cross(e2,v3_t);
          v2_t=cross(v3_t,v1)+cross(v3,v1_t);
       end
      

       
     
     
