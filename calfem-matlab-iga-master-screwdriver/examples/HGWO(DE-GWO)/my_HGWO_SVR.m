%  function [bestc,bestg,test_pre]=my_HGWO_SVR(para,input_train,output_train,input_test,output_test)
% 参数向量 parameters [n,N_iteration,beta_min,beta_max,pCR]
% n为种群规模，N_iteration为迭代次数
% beta_min 缩放因子下界 Lower Bound of Scaling Factor
% beta_max=0.8; % 缩放因子上界 Upper Bound of Scaling Factor
% pCR 交叉概率 Crossover Probability
% 要求输入数据为列向量（矩阵）
% tic


load initial 

% [INBz1,INNz1,Br1,Bz11,Bz] =Midgridx(B,lb,ub,Xmax,Xmin,KV,deg) ;
para=[50,600,0.2,0.8,0.2];
%% 数据归一化
% [input_train,rule1]=mapminmax(input_train');
% [output_train,rule2]=mapminmax(output_train');
% input_test=mapminmax('apply',input_test',rule1);
% output_test=mapminmax('apply',output_test',rule2);
% input_train=input_train';
% output_train=output_train';
% input_test=input_test';
% output_test=output_test';
bB=cell2mat(B1);
B01=B1;
 input_weight=bB(4,:,:);
% [input_weight,rule1]=mapminmax(input_weight');
% [input_train,rule1]=mapminmax(input_DATA');

%% 利用差分进化-灰狼优化混合算法（DE_GWO）选择最佳的SVR参数
nPop=para(1); % 种群规模 Population Size
MaxIt=para(2); % 最大迭代次数Maximum Number of Iterations
%  nVar=2; % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
 nVar=(size(B1,2))*size(B1,3); % 自变量维数，此例需要优化两个参数c和g Number of Decision Variables
VarSize=[1,nVar]; % 决策变量矩阵大小 Decision Variables Matrix Size
beta_min=para(3); % 缩放因子下界 Lower Bound of Scaling Factor
beta_max=para(4); % 缩放因子上界 Upper Bound of Scaling Factor
pCR=para(5); %  交叉概率 Crossover Probability
% lb=[0.01,0.01]; % 参数取值下界
% ub=[100,100]; % 参数取值上界
lb=0.8*ones(1,nVar); % 参数取值下界
ub=1.2*ones(1,nVar);% 参数取值上界                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
%% 初始化
% 父代种群初始化
tic
[INBz1,INNz1,Br1,Bz11,Bz] =Midgridx(B0,lb,ub,Max,Min,KV,deg) ;
[ Pp1,nsub,ssub,idsub]=suoyin_1(INBz1,INNz1,Br1,Bz11,B01,Bz); %1初始,2参数化修改
INPUT_=[];OUTPUT_=[];
parent_Position=[1.10974134805838,0.955472553599478,1.09064566978359,1.04595113324898,0.824072623171227,0.895852350299567,0.811204483909518,1.01413346635150,1.00183675614083,1.11437434727461,1.08802616224102,0.980070185339065,0.829262855188528,0.897671367848784,1.06908333609771,1.12661114305046,0.994187685386687,1.19918508174534,0.923973025628442,0.966018939637052,1.09281385892947,0.855840421356620,1.12868045779865,0.945091179903625,0.976780043351720,0.827494919970892,0.835245858172596,0.855405267875801,0.926841949990231,1.05724201221350,0.929585966182736,0.859299958128166,0.820748749192408,0.924436501358333,1.08712010001684,1.10467643865676,1.04983567440254,0.962737490310485,1.02612321489979,1.10951115389376,0.852606078229835,0.960615286558764,1.08706564020660,0.919229058561365,0.893246560777172,1.14136976286664,1.17912355654492,0.980961954087304,1.05993393777085,0.951959043181385,0.895502429260318,1.02834763725118,0.982007867748141,0.964882748393323,1.09947636684100,1.09191618435597,0.980841927303585,1.00031262369919,1.01828366453379,1.12654581607200,1.11252536602078,1.14671875021614,1.10425187821848,1.16890306219480,1.01332895010013,1.17668140373899,0.929526859910169,1.18293734678843,1.07780077287652,1.14761884400653,1.04252687981535,1.15963618592642,1.10629551467828,1.01693018872767,0.893963467011676,0.981667742810573,1.15688765444471,0.899661387759560,0.983687454167379,1.14647229594946,1.15877251858532,0.871805545305552,0.945548679361160,1.08361923693259,0.889568234791084,1.05410329723913,0.894720038796099,1.08371536217753,0.800192843943794,0.913358325537637,0.853424566925211,1.13420661689393,1.16833281609482,1.14817849114297,1.02929346156815,0.955693612680740,0.849146495801613,1.16316675681106,1.16058846401744,0.873619551367493,0.892328381550288,1.03808309912340,1.11963165631083,1.16600572604383,0.970553170594207,0.833972731758118,0.843954173988216,1.13809732486006,0.950272728854258,1.13744092169703,0.967482155808490,1.04373212160118,1.10946929630944,1.12540644497764,0.944607996848662,1.06277587961074,1.09018082020271,0.830701176992968,1.00364198982298,0.834177628009994,0.840722769380936,0.829774529065487,0.866527094769439,1.08453240593539,0.999937619205205,1.16081767379792,1.00603607649835,0.932558679878727,1.07758197404274,0.823150721606955,1.00098635123565,0.844138238673563,1.04557869386991,1.14954170700145,0.963942868455046,0.959169332030724,0.918568540652835,1.09334579705219,1.08934370201609,0.971853829997869,0.973739755331850,0.934086425713125,1.15330379170754,1.19191309625087,0.802037827969611,1.19144143735636,1.15872152091188,1.05596914889594,0.840273155707395,1.13834784258034,0.934039892177727,1.06081469918108,0.900344178433272,1.01159059956636,0.914917114049558,0.936089780038121,0.985175685972428,0.937783620428540,1.10785885526163,1.03660806237531,1.02628256154821,0.838593085426555,1.03275890052339,0.944185060653206,0.909514976859739,0.949936066410436,1.17286392210351,0.981276376820678,0.838408209342747,0.920270435961839,1.19201763177805,0.992834195340651,1.16176394815700,0.925893495101783,1.12403901036766,1.08000483701579,0.806194525574507,1.00919744817415,1.09627085151554,1.12804278397840];
 parent_Val= objecfun(parent_Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
low_efficiency=xlsread('grade2.xls');
 Rank=randi(size(low_efficiency,2),1,size(low_efficiency,2));
 matlabpool local 4;  
for il=1:numel(Rank)/4
    il
    
    as=size(low_efficiency,1);
   
    lb_n=0.95*low_efficiency(2:as,il)';
    lb_n(lb_n<0.8)=0.8;
    ub_n=1.05*low_efficiency(2:as,il)';
    ub_n(ub_n>1.2)=1.2;
    
 inter=10;
 parent_Val=zeros(nPop,inter); % 目标函数值
 
 %% 并行随机
  
for js=1:inter
    js
  parent_Position=init_individual(lb_n,ub_n,nVar,nPop); % 随机初始化位置
% parent_Position=0.5*ones(1,nVar);
 for is=1:nPop % 遍历每个个体
% for i=1:1 % 遍历每个个体
%     parent_Val(i)=fobj(parent_Position(i,:),input_train,output_train,input_test,output_test); % 计算个体目标函数值
    parent_Val(is,js)= objecfun(parent_Position(is,:),nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
 end
 input{js}=parent_Position;
end

INPUT_NORML=reshape((cell2mat(input))',nVar,numel(cell2mat(input))/nVar);
OUT_NORMAL=reshape(parent_Val',1,numel(parent_Val));
INPUT_=[INPUT_; INPUT_NORML'];
OUTPUT_=[OUTPUT_; OUT_NORMAL'];
end
matlabpool close
toc
%{
% Rank=randi(size(INPUT_NORMAL,1),1,size(INPUT_NORMAL,1));
in_low=INPUT_NORMAL;
out_low=OUTPUT_NORMAL;
%反向学习方法对初始种群取优
parent_Position_fan=(lb(1)+ub(1)).*ones(size(parent_Position))-parent_Position;
for is=1:nPop % 遍历每个个体取反
%     parent_Val(i)=fobj(parent_Position(i,:),input_train,output_train,input_test,output_test); % 计算个体目标函数值
    parent_Val_fan(is)= objecfun(parent_Position_fan(is,:),nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
end
parent_Val=[parent_Val ; parent_Val_fan'];
parent_Position=[parent_Position; parent_Position_fan];

% 突变种群初始化
mutant_Position=init_individual(lb,ub,nVar,nPop); % 随机初始化位置
mutant_Val=zeros(nPop,1); % 目标函数值
for is=1:nPop % 遍历每个个体
%     mutant_Val(i)=fobj(mutant_Position(i,:),input_train,output_train,input_test,output_test); % 计算个体目标函数值
     mutant_Val(is)=objecfun(mutant_Position(is,:),nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1); 
end
% 子代种群初始化
child_Position=init_individual(lb,ub,nVar,nPop); % 随机初始化位置
child_Val=zeros(nPop,1); % 目标函数值
for is=1:nPop % 遍历每个个体
%     child_Val(i)=fobj(child_Position(i,:),input_train,output_train,input_test,output_test); % 计算个体目标函数值
    child_Val(is)=objecfun(child_Position(is,:),nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1); 
end

mutant_Val=zeros(nPop,1); % 目标函数值
child_Val=zeros(nPop,1); % 目标函数值
%% 确定父代种群中的Alpha,Beta,Delta狼
[~,sort_index]=sort(parent_Val); % 父代种群目标函数值排序
parent_Alpha_Position=parent_Position(sort_index(1),:); % 确定父代Alpha狼
parent_Alpha_Val=parent_Val(sort_index(1)); % 父代Alpha狼目标函数值
parent_Beta_Position=parent_Position(sort_index(2),:); % 确定父代Beta狼
parent_Delta_Position=parent_Position(sort_index(3),:); % 确定父代Delta狼
parent_Val=parent_Val(sort_index(1:numel(sort_index)/2));
parent_Position=parent_Position(sort_index(1:numel(sort_index)/2),:);
[~,sort_index]=sort(parent_Val);
%% 迭代开始
BestCost=zeros(1,MaxIt);
BestCost(1)=parent_Alpha_Val;
% figure
for it=1:MaxIt
    a=2-it*((2)/MaxIt); % 对每一次迭代，计算相应的a值，a decreases linearly fron 2 to 0
    a1=2*exp(((it/MaxIt)^2)*log(0.01/2));  %a=amax*exp(((it/MaxIt)^eta_alpha)*ln(amin/amax)),
    a3=2*exp(((it/MaxIt)^3)*log(0.01/2));   %a=amax*exp(((it/MaxIt)^eta_egma)*ln(amin/amax)),
    a2=(a1+a3)/2;
%     if it>=2*MaxIt/3
    % 更新父代个体位置
    parent_Val1=parent_Val;
    parent_Position1=parent_Position;
    for par=1:nPop % 遍历父代个体

%       if rand(0,1)<0.7  %趋优算子设置为0.5
        for var=1:nVar % 遍历每个维度            
            % Alpha狼Hunting
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]            
%             A1=2*a*r1-a; % 计算系数A
            A1=2*a1*r1-a1; % 计算系数A
            C1=2*r2; % 计算系数C
            D_alpha=abs(C1*parent_Alpha_Position(var)-parent_Position1(par,var));
            X1=parent_Alpha_Position(var)-A1*D_alpha;
            % Beta狼Hunting
            r1=rand();
            r2=rand();            
%             A2=2*a*r1-a; % 计算系数A
            A2=2*a2*r1-a2; % 计算系数A
            C2=2*r2; % 计算系数C
            D_beta=abs(C2*parent_Beta_Position(var)-parent_Position1(par,var));
            X2=parent_Beta_Position(var)-A2*D_beta;
            % Delta狼Hunting
            r1=rand();
            r2=rand();
%             A3=2*a*r1-a; % 计算系数A
            A3=2*a3*r1-a3; % 计算系数A
            C3=2*r2; % 计算系数C
            D_delta=abs(C3*parent_Delta_Position(var)-parent_Position1(par,var));
            X3=parent_Delta_Position(var)-A3*D_delta;
            % 位置更新，防止越界
%              X=(X1+X2+X3)/3;
            %利用权重系数更新
            wnum=sum(parent_Val1(sort_index(1:3)));
            if parent_Val1(par)>=wnum/3
            X=parent_Val1(sort_index(1))/wnum*X1+parent_Val1(sort_index(2))/wnum*X2+parent_Val1(sort_index(3))/wnum*X3;
            X=0.5*X+0.5*rand*(parent_Delta_Position(var)-parent_Position1(par,var));  %考虑个体记忆功能
            
            else
            X=(X1+X2+X3)/3;  
            end
             
            X=max(X,lb(var));
            X=min(X,ub(var));
            parent_Position(par,var)=X;
 
        end
%       else
%           for var=1:nVar % 遍历每个维度    
%           A=randperm(nVar); % 个体顺序重新随机排列
%           a=A(1);
%             X=(parent_Alpha_Position(a)+parent_Beta_Position(a)+parent_Delta_Position(a))/3;
%             X=max(X,lb(var));
%             X=min(X,ub(var));
%             parent_Position(par,var)=X;
%           end
%         parent_Val(par)=fobj(parent_Position(par,:),input_train,output_train,input_test,output_test); % 计算个体目标函数值
%       end
       parent_Val(par)= objecfun(parent_Position(par,:),nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1); 
       
       %引入个体记忆后对历史最优进行比较
       if parent_Val(par)> parent_Val1(par)
           parent_Position(par,:)=parent_Position1(par,:);
           parent_Val(par)=parent_Val1(par);
       end
    end
%     else
    % 产生变异（中间体）种群
    for mut=1:nPop
        A=randperm(nPop); % 个体顺序重新随机排列
        A(A==mut)=[]; % 当前个体所排位置腾空（产生变异中间体时当前个体不参与）
        a=A(1);
          b=A(2);
          c=A(3);
        beta=unifrnd(beta_min,beta_max,VarSize); % 随机产生缩放因子
        y=parent_Position(a)+beta.*(parent_Position(b)-parent_Position(c)); % 产生中间体
        % 防止中间体越界
        y=max(y,lb);
		y=min(y,ub);
        mutant_Position(mut,:)=y;
    end
    % 产生子代种群，交叉操作 Crossover
    for child=1:nPop
        x=parent_Position(child,:);
        y=mutant_Position(child,:);
        z=zeros(size(x)); % 初始化一个新个体
        j0=randi([1,numel(x)]); % 产生一个伪随机数，即选取待交换维度编号？？？
        for var=1:numel(x) % 遍历每个维度
            if var==j0 || rand<=pCR % 如果当前维度是待交换维度或者随机概率小于交叉概率
                z(var)=y(var); % 新个体当前维度值等于中间体对应维度值
            else
                z(var)=x(var); % 新个体当前维度值等于当前个体对应维度值
            end
        end
        child_Position(child,:)=z; % 交叉操作之后得到新个体
%         child_Val(child)=fobj(z,input_train,output_train,input_test,output_test); % 新个体目标函数值
        child_Val(child)= objecfun(z,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1); 
    end
    % 父代种群更新
    for par=1:nPop
        if child_Val(par)<parent_Val(par) % 如果子代个体优于父代个体
            parent_Val(par)=child_Val(par);    %更新父代个体特征值
            parent_Position(par,:)=child_Position(par,:);  % 更新父代个体
        end
    end
%     end
    % 确定父代种群中的Alpha,Beta,Delta狼
    [~,sort_index]=sort(parent_Val); % 父代种群目标函数值排序
    parent_Alpha_Position=parent_Position(sort_index(1),:); % 确定父代Alpha狼
    parent_Alpha_Val=parent_Val(sort_index(1)); % 父代Alpha狼目标函数值
    parent_Beta_Position=parent_Position(sort_index(2),:); % 确定父代Beta狼
    parent_Delta_Position=parent_Position(sort_index(3),:); % 确定父代Delta狼
    BestCost(it)=parent_Alpha_Val;
     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
%     PlotCosts1(BestCost,it);
%     pause(0.01);
%     hold on
end
bestc=parent_Alpha_Position(1,1);
bestg=parent_Alpha_Position(1,2);
%% 图示寻优过程
figure ( )
plot(BestCost);
xlabel('Iteration');
ylabel('Best Val');
grid on;

%{
%% 利用回归预测分析最佳的参数进行SVM网络训练
cmd_cs_svr=['-s 3 -t 2',' -c ',num2str(bestc),' -g ',num2str(bestg)];
model_cs_svr=svmtrain(output_train,input_train,cmd_cs_svr); % SVM模型训练
%% SVM网络回归预测
[output_test_pre,~]=svmpredict(output_test,input_test,model_cs_svr); % SVM模型预测及其精度
test_pre=mapminmax('reverse',output_test_pre',rule2);
test_pre = test_pre';
%}
%}