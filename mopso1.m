clc;
clear;
close all;
format long
%% 
global  nu0 t0 tf E_int E_eq g k tau Ein Ef
%%
% E_int=input('plz enter Eint_exp');
% E_eq=input('plz enter Eeq_exp');
% nu0=input('plz enter Poisson’s ratio');
% disp('----------------------------------')
%%
E_int=5.09e6;
E_eq=0.78e6;
nu0=0.3;
t0=0;
tf=1000;
T=[t0,tf];
syms t
%%
CostFunction=@(x) GK(x);      
%
nVar=3;            
%
VarSize=[1 nVar];    
%
VarMin=1;            
VarMax=100;     
%%
MaxIt=5;          
%
nPop=50;           
%
nRep=50;           
%% 
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt((phi^2)-4*phi));
w=chi;
wdamp=1;
c1=chi*phi1;
c2=chi*phi2;
%% 
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
%
nGrid=20;          
alpha=0.1;          
%
beta=.5;            
gamma=2;          
%
mu=0.1;             
%% 
empty_particle.Position=[];
empty_particle.Velocity=[];
%
empty_particle.g=[];
empty_particle.k=[];
empty_particle.tau=[];
empty_particle.Eopt_in=[];
empty_particle.Eopt_eq=[];
%
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
%
empty_particle.Best.g=[];
empty_particle.Best.k=[];
empty_particle.Best.tau=[];
empty_particle.Best.Eopt_in=[];
empty_particle.Best.Eopt_eq=[];
%
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];
pop=repmat(empty_particle,nPop,1);
%%
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Velocity=zeros(VarSize);
    % 
    pop(i).g=g;
    pop(i).k=k;
    pop(i).tau=tau;
    pop(i).Ein=Ein;
    pop(i).Ef=Ef;
    %
Cost1=CostFunction(pop(i).Position);
pop(i).Cost=Cost1;
   %
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
    pop(i).Best.g=g;
    pop(i).Best.k=k;
    pop(i).Best.tau=tau;
    pop(i).Best.Ein=Ein;
    pop(i).Best.Ef=Ef;
%                 
end
% 
pop=DetermineDomination(pop);
rep=pop(~[pop.IsDominated]);
Grid=CreateGrid(rep,nGrid,alpha);
for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end
%%
for it=1:MaxIt
    for i=1:nPop
        leader=SelectLeader(rep,beta);
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        %% 
        pop(i).Velocity=...
            max(pop(i).Velocity,VelMin);
        pop(i).Velocity=...
            max(pop(i).Velocity,VelMax);
        %% 
        IsOutside=...
            (pop(i).Position<VarMin | pop(i).Position>VarMax);
        pop(i).Velocity(IsOutside)=...
            -pop(i).Velocity(IsOutside);
    %% 
        pop(i).Position=...
            max(pop(i).Position,VarMin);
        pop(i).Position=...
            min(pop(i).Position,VarMax);
        %
        Cost2 = CostFunction(pop(i).Position);
        pop(i).Cost =Cost2;
        % 
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
        Cost3=CostFunction(NewSol.Position);
        NewSol.Cost=Cost3;
       %
        if Dominates(NewSol,pop(i))
            pop(i).Position=NewSol.Position;
            pop(i).Cost=NewSol.Cost; 
          %  
        elseif Dominates(pop(i),NewSol)
            % 
        else
            if rand<0.5
                pop(i).Position=NewSol.Position;
                pop(i).Cost=NewSol.Cost;
            end
        end
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
%  
            pop(i).Best.g=pop(i).g;
            pop(i).Best.k=pop(i).k;
            pop(i).Best.tau=pop(i).tau;
            pop(i).Best.Eopt_in=Ein;  
            pop(i).Best.Eopt_eq=Ef;  
        elseif Dominates(pop(i).Best,pop(i))
            % 
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
              %
                pop(i).Best.g=pop(i).g;
                pop(i).Best.k=pop(i).k;
                pop(i).Best.tau=pop(i).tau;
                pop(i).Best.Eopt_in=Ein;
                pop(i).Best.Eopt_eq=Ef; 
                end
            end
        end
    rep=[rep;pop(~[pop.IsDominated])];
    % 
    rep=DetermineDomination(rep);
    %
    rep=rep(~[rep.IsDominated]);
    %
    Grid=CreateGrid(rep,nGrid,alpha);
%
%%
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
% 
    end
% 
    if numel(rep)>nRep
        %
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        %
    end
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    % 
    w=w*wdamp;
    end   
%% 
disp('--------------------------------')
    disp(['g = ',num2str(rep(end).Best.g)])
    disp(['k = ',num2str(rep(end).Best.k)])
    disp(['tau = ',num2str(rep(end).Best.tau)])
    disp(['Eoptimal_init (Pa)= ',num2str(rep(end).Best.Eopt_in)])
    disp(['Eoptimal_equal (Pa)= ',num2str(rep(end).Best.Eopt_eq)])