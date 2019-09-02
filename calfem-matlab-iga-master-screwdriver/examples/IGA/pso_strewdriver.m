clc;
clear;
close all;
tic
load initial 
bB=cell2mat(B1);
B01=B1;
 %{
% for tubular structure
 nVar=(size(B1,2)-2)*size(B1,3);
lb=0.8*ones(1,nVar); % 参数取值下界
ub=1.2*ones(1,nVar);% 参数取值上界  
%}

nVar1 = 1;    % Number of Decision Variables
nVar2 = 6;
nVar3 = 3;
nVar=nVar1+nVar2+nVar3;
VarSize1 = [1 nVar1];
VarSize2 = [1 nVar2];
VarSize3 = [1 nVar3];


% %% Problem Definition
% 
% CostFunction=@(x) REiga_pro3(x);        % Cost Function
% 
% % nVar=1;            % Number of Decision Variables
% % 
% % VarSize=[1 nVar];   % Size of Decision Variables Matrix
% 
% nVar1 = 1;    % Number of Decision Variables
% 
% % nVar2 = 1;
%  VarSize2 = [1 nVar1];
% VarSize1 = [1 nVar1];
% 
% % VarMin=1.45;        % Lower Bound of Variables
% % VarMax=1.8;         % Upper Bound of Variables
% 
V1Min = 0;   % Lower Bound of Variables
V1Max = 2;    % Upper Bound of Variables
V2Min = 0.8;   % Lower Bound of Variables
V2Max = 1.2;
V3Min = 0.5;   % Lower Bound of Variables
V3Max = 1.5;
VarMax=[V1Max.*ones(VarSize1),V2Max.*ones(VarSize2),V3Max.*ones(VarSize3)];
VarMin=[V1Min.*ones(VarSize1),V2Min.*ones(VarSize2),V3Min.*ones(VarSize3)];

ub=VarMax;
lb=VarMin;

[INBz1,INNz1,Br1,Bz11,Bz] =Midgridx2(B0,lb,ub,KV0,deg) ;    %Midgridx2 for screwdriver
[ Pp1,nsub,ssub,idsub]=suoyin_1(INBz1,INNz1,Br1,Bz11,B01,Bz); %1初始,2参数化修改
% 
% KV0
%% PSO Parameters
% load 93
% GlobalBest=PSO{numel(PSO)};
% x=GlobalBest(1:(numel(GlobalBest)-1));
% tic
% [objtive,stress_n,aT]=obj_strewdriver(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
% toc

MaxIt=100;      % Maximum Number of Iterations

nPop=40;        % Population Size (Swarm Size)

% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
% VelMax=0.1.*(ub-lb);
% VelMin=-VelMax;
VelMax1=0.1*(V1Max-V1Min);
VelMin1=-VelMax1;
VelMax2=0.1*(V2Max-V2Min);
VelMin2=-VelMax2;
VelMax3=0.1*(V3Max-V3Min);
VelMin3=-VelMax3;
VelMin=[VelMin1.*ones(VarSize1),VelMin2.*ones(VarSize2),VelMax3.*ones(VarSize3)];
VelMax=[VelMax1.*ones(VarSize1),VelMax2.*ones(VarSize2),VelMin3.*ones(VarSize3)];
%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
%     x=init_individual(lb,ub,nVar,nPop); % 随机初始化位
    x1=unifrnd(V1Min,V1Max,VarSize1);
    x2=unifrnd(V2Min,V2Max,VarSize2);
    x3=unifrnd(V3Min,V3Max,VarSize3);
% xRange=ub-lb;
% xLower=lb;
% x=rand(1,nVar).*xRange+xLower;
%     
     x=[x1,x2,x3];
    particle(i).Position=x;
    % Initialize Velocity
    v1=zeros(VarSize1);v2=zeros(VarSize2);v3=zeros(VarSize3);v=[v1,v2,v3];
%     v=zeros(1,nVar);
    particle(i).Velocity=v;
    % Evaluation
%     particle(i).Cost=CostFunction(particle(i).Position);
 
% profile on 
   particle(i).Cost=obj_strewdriver(particle(i).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
%     profile viewer
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(1,nVar).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(1,nVar).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
%         particle(i).Cost = CostFunction(particle(i).Position);
        particle(i).Cost =obj_strewdriver(particle(i).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
        
    end
    
    BestCost(it)=GlobalBest.Cost;
%     PSO{it}(1:nVar)=GlobalBest.Position;
%     PSO{it}(nVar+1)=GlobalBest.Cost;
  
    bestposition_pso(it,:)=GlobalBest.Position;
    bestcost_pso=BestCost(it);
    PSO_ERROR(it)=error_objecfun(GlobalBest.Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
%         % Plot F1 Costs
%     figure(1);
%     PlotCosts(GlobalBest,it);
%     pause(0.01);
%     hold on
    w=w*wdamp;
    save it
    name1='it.mat';
    na='.mat';
    name2= strcat(num2str(it),na);
    movefile(name1,name2);
end
toc
BestSol = GlobalBest;

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

