%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

%


clc;
clear;
close all;

%% Problem Definition

CostFunction=@(x) REiga_pro3(x);    % Cost Function

nVar=1;            % Number of Decision Variables

VarSize1=[1 nVar];   % Decision Variables Matrix Size
VarSize2=[1 nVar];   % Decision Variables Matrix Size
% VarMin=1.15;          % Lower Bound of Decision Variables  %BENCHMARK
% VarMax=1.4;          % Upper Bound of Decision Variables
VarMin1=1.2;          % Lower Bound of Decision Variables  %SCREWDRIVER
VarMax1=1.35;          % Upper Bound of Decision Variables
VarMin2=1.9;          % Lower Bound of Decision Variables  %SCREWDRIVER
VarMax2=2.2;          % Upper Bound of Decision Variables


%% DE Parameters

MaxIt=40;      % Maximum Number of Iterations

nPop=40;        % Population Size

beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor

pCR=0.2;        % Crossover Probability

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop

    x1=unifrnd(VarMin1,VarMax1,VarSize1);
    x2=unifrnd(VarMin2,VarMax2,VarSize1);
    x=[x1,x2];
    pop(i).Position=x;
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        beta1=unifrnd(beta_min,beta_max,VarSize1);
        beta1=beta1';
        y1=pop(a).Position+beta1.*(pop(b).Position-pop(c).Position);
        y1 = max(y1, VarMin1);     %judging condition
		y1 = min(y1, VarMax1);
        beta2=unifrnd(beta_min,beta_max,VarSize1);
        beta2=beta2';
        y2=pop(a).Position+beta2.*(pop(b).Position-pop(c).Position);
        y2 = max(y2, VarMin2);     %judging condition
		y2 = min(y2, VarMax2);
        y=[y1,y2];
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCostde(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCostde(it))]);
    
end

%% Show Results

figure;
%plot(BestCost);
semilogy(BestCostde, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
