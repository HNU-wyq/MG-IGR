
function Kee =unit_assemble_stiff(Km,Kb,nen)

for m = 1:nen
    for n = 1:nen                   %局部坐标下的平面单元刚度矩阵与弯曲刚度矩阵共同组装刚度矩阵
        Kee((m-1)*6+1:(m-1)*6+3,(n-1)*6+1:(n-1)*6+3) = Km((m-1)*3+1:(m-1)*3+3,(n-1)*3+1:(n-1)*3+3); %平面刚度矩阵
        Kee((m-1)*6+4:(m-1)*6+6,(n-1)*6+4:(n-1)*6+6) = Kb((m-1)*3+1:(m-1)*3+3,(n-1)*3+1:(n-1)*3+3);  %弯曲刚度矩阵
        %Kee((i-1)*6+6,(j-1)*6+6) = 0;       %ke（6i,6j）元素全为0   

    end
end