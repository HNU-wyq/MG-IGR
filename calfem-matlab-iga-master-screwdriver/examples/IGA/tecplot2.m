function  tecplot2(IEN,gcoord,stress_n,a,Vz)

nel = size(IEN,2);
node= size(gcoord,1);
nodedata=fopen('Mise.plt','w');               
fprintf(nodedata,'TITLE="data"\n');

fprintf(nodedata,'VARIABLES=,"X", "Y","Z" ,"sig1","sig2","sig3","ux","uy","uz","qz"\n');

fprintf(nodedata,'ZONE T="%d  "  ,  N=%d, E=%d, ET=Quadrilateral, F=FEPOINT\n',1,node,nel); % 注意这里的双引号
sigm1=stress_n(:,1);sigm2=stress_n(:,2);sigm3=stress_n(:,3);
ux=a(1:6:size(a));uy=a(2:6:size(a));uz=a(3:6:size(a)); qz=Vz(:,3);

for i=1:node
    fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n',...
      gcoord(i,1),gcoord(i,2),gcoord(i,3),sigm1(i),sigm2(i),sigm3(i),ux(i),uy(i),uz(i),qz(i));
%      fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f\n',gcoord(i,1),gcoord(i,2),gcoord(i,3),tempfull(i));
end  

for j=1:nel
    for i=1:size(IEN,1)
        fprintf(nodedata,'%d       ',IEN(i,j));  
    end
    fprintf(nodedata,'\n');
end
end