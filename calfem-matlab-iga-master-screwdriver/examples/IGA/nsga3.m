% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

clc;
clear all;
close all;

%% Problem Definition
tic

load initial 
bB=cell2mat(B1);
B01=B1;

nVar1 = 1;    % Number of Decision Variables
nVar2 = 6;
nVar3 = 3;
nVar=nVar1+nVar2+nVar3;
VarSize1 = [1 nVar1];
VarSize2 = [1 nVar2];
VarSize3 = [1 nVar3];

V1Min = 0;   % Lower Bound of Variables
V1Max = 2;    % Upper Bound of Variables
V2Min = 0.8;   % Lower Bound of Variables
V2Max = 1.2;
V3Min = 0.5;   % Lower Bound of Variables
V3Max = 1.5;
x1=unifrnd(V1Min, V1Max, VarSize1);
x2=unifrnd(V2Min, V2Max, VarSize2);
x3=unifrnd(V3Min, V3Max, VarSize3);
x=[x1,x2,x3];

VarMax=[V1Max.*ones(VarSize1),V2Max.*ones(VarSize2),V3Max.*ones(VarSize3)];
VarMin=[V1Min.*ones(VarSize1),V2Min.*ones(VarSize2),V3Min.*ones(VarSize3)];

ub=VarMax;
lb=VarMin;

[INBz1,INNz1,Br1,Bz11,Bz] =Midgridx2(B0,lb,ub,KV0,deg) ;    %Midgridx2 for screwdriver
[ Pp1,nsub,ssub,idsub]=suoyin_1(INBz1,INNz1,Br1,Bz11,B01,Bz); %1初始,2参数化修改


CostFunction = @(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1) obj_strewdriver(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);  % Cost Function
% nVar = 1;    % Number of Decision Variables
%  nVar2 = 1;
%  nVar3 = 2;

% VarSize1 = [1 nVar]; % Size of Decision Variables Matrix
% VarSize2 = [1 nVar];
% % % VarSize3 = [1 nVar3];
% 
% % V1Min = 1.15;   % Lower Bound of Variables    %benchmark
% % V1Max = 1.4;    % Upper Bound of Variables
% V1Min = 1.18;   % Lower Bound of Variables    %screwdriver
% V1Max = 1.26;    % Upper Bound of Variables
% V2Min = 1.9;   % Lower Bound of Variables    %screwdriver
% V2Max = 2.3;    % Upper Bound of Variables

% x1=unifrnd(V1Min, V1Max, VarSize1);
% x2=unifrnd(V2Min, V2Max, VarSize2);
% x=[x1,x2];
% % Number of Objective Functions
% %nObj = numel(CostFunction(unifrnd(V1Min, V1Max, VarSize1)));  %ladingcouyang
 nObj = numel(CostFunction(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1));



%% NSGA-II Parameters

% Generating Reference Points
nDivision = 10;
Zr = GenerateReferencePoints(nObj, nDivision);

MaxIt = 100;  % Maximum Number of Iterations 100

nPop = 40;  % Population Size  40

pCrossover = 0.8;       % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2); % Number of Parnets (Offsprings)

pMutation = 0.3;       % Mutation Percentage
nMutation = round(pMutation*nPop);  % Number of Mutants

mu = 0.02;     % Mutation Rate

%sigma = 0.1*(VarMax-VarMin); % Mutation Step Size
sigma = 0.1*(V1Max-V1Min);

%% select Parameters

params.nPop = nPop;
params.Zr = Zr;
params.nZr = size(Zr,2);
params.zmin = [];
params.zmax = [];
params.smin = [];

%% Initialization

disp('Staring NSGA-III ...');

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.NormalizedCost = [];
empty_individual.AssociatedRef = [];
empty_individual.DistanceToAssociatedRef = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    %pop(i).Position = unifrnd(VarMin, VarMax, VarSize); % varibles
     x1=unifrnd(V1Min, V1Max, VarSize1);
     x2=unifrnd(V2Min, V2Max, VarSize2);
     x3=unifrnd(V3Min, V3Max, VarSize3);
     x=[x1 x2 x3];
%  x=x1;
    pop(i).Position = x;
    pop(i).Cost = CostFunction(pop(i).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
end
clear x1 x2 x3
% Sort Population and Perform Selection
[pop, F, params] = SortAndSelectPopulation(pop, params);

disp('Initialization.');
%% NSGA-II Main Loop

for it = 1:MaxIt
 
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2
        k
        i1 = randi([1 nPop]);
        p1 = pop(i1);

        i2 = randi([1 nPop]);
        p2 = pop(i2);

        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position,VarMax,VarMin);

        popc(k, 1).Cost = CostFunction(popc(k, 1).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);

    end
    popc = popc(:);

    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation

        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, mu, sigma,VarMax,VarMin);

        popm(k).Cost = CostFunction(popm(k).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);

    end
 
    % Merge
    pop = [pop
           popc
           popm]; %#ok
    
    % Sort Population and Perform Selection
    [pop, F, params] = SortAndSelectPopulation(pop, params);
    
    % Store F1
%     F1 = pop(F{1});
    if numel(F{1})==1
    F1 = pop(F{1});
    else
     F1 = pop(F{1}(1)); 
    end
    
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(F1.Cost)]);

    % Plot F1 Costs
    figure(1);
    PlotCosts(F1,it);
    pause(0.01);
    hold on
    GA(it,1:nVar)=F1.Position;
    GA(it,nVar+1)=F1.Cost;
    GA_Error(it)=error_objecfun(F1.Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
end
% reana(F1);

%% Results

disp(['Final Iteration: Number of F1 Members = ' num2str(numel(F1))]);
disp('Optimization Terminated.');
toc

