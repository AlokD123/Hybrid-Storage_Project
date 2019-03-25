%Storage sizing test with ALP optimization
%---------- Using APPROX LP ONLY!! -----------
%ADDED Regenerative Braking
%Realistic sizing

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[4,4];
min_E_SIZE=[1,1];

%Parameters
SCALE_BATT=89/1000*379/max_E_SIZE(1); %kWh/kg*kg/gridpt for battery
SCALE_SC=4/1000*4/max_E_SIZE(2); %kWh/kg*kg/gridpt for supercap

INFCOST=1e6;

global E_MAX; global E_MIN;


size_iter=0; %Storage size iteration
PF_opt_mtx=[]; %Transition weights mtx
%optNextE_arr=[]; %Next states array
g_opt_vect=[]; %Stage costs vector
%Store transition weights and stage costs for optimal policies for EACH iteration
PF_opt={}; g_opt={}; Exp_CostToGo={}; optCost_size={}; g_opt_mtx={}; Exp_CostToGo_mtx={};

%Optimal size ratio and costs for each value of C1,C2
optRatio=[]; optCosts={};

global RES_E1; global RES_E2; global RES_L; global RES_U1; 

%Step sizes
maxE_stepSize_E1=1;
maxE_stepSize_E2=1;

%CONSTANT resolutions in simulation
RES_E1=1/(7500*2); %*maxE_stepSize_E1);%*size1_mult);
RES_E2=1/(330*2); %(maxE_stepSize_E2);%size2_mult);
RES_L=1/4; %4/(2*size1_mult+30*size1_mult);
RES_U1=1; %/size1_mult;

mult_cost_idx=0;

for cost_mult_1=1:5
    for cost_mult_2=1:5    
        
    mult_cost_idx=mult_cost_idx+1;
               
    %Cost for size
    c1_fact=RES_E1*10; c1=c1_fact*(cost_mult_1-1); %0.01;
    c2_fact=RES_E2/20; c2=c2_fact*(cost_mult_2-1); %0.5;
    
    %(Re-)initialize
    minCost=Inf; %Default: initial value for J(initial state)+cS cost
    vectS_netOptVal=[];%Store CONSTANT state optimal value for each given S
    totCost=[]; %Store total cost for each size (convex function)
    
    %Counters for grid
    E1_counter=length(min_E_SIZE(1):maxE_stepSize_E1:max_E_SIZE(1))+1;
    for max_E1=fliplr(min_E_SIZE(1):maxE_stepSize_E1:max_E_SIZE(1))
        E1_counter=E1_counter-1;

        E2_counter=length(min_E_SIZE(1):maxE_stepSize_E2:max_E_SIZE(2))+1;
        for max_E2=fliplr(min_E_SIZE(1):maxE_stepSize_E2:max_E_SIZE(2))
            E2_counter=E2_counter-1;

            %Go through feasible set for E_SIZE
            size1_mult=max_E1; size2_mult=max_E2;

            E_MAX=[7500*2*size1_mult;330*2*size2_mult];

            size_iter=size_iter+1; %Next size up

            %if size_iter<95
        
            %else
            
            ApproxLP_sol_IHDP_v18; %Get optimal values for this size
            %GetCtrlPolicy_OptQVals_v2; %Get optimal policy matrix

            %Store optimal values, for reference;
            optVal_size{size_iter}=ConvCosts; %<--------- DIFFERENCE IS LARGE

            %Get optimal value for MINIMUM allowable capacity w/ ZERO LOAD.
            %NOTE: have index 2 ASSUMING E_MIN=[0,0] is a state.
            optVal_initE=ConvCosts((min_E_SIZE(1)-0)+1,(min_E_SIZE(1)-0)+1,round((-MIN_LOAD)*RES_L+1));

            vectS_netOptVal=[vectS_netOptVal;optVal_initE]; %Store in vector

            %Get optimal storage size till this point
            if (optVal_initE + c1*E_MAX(1) + c2*E_MAX(2)) < minCost   %c is COST-PER-UNIT, **not** cost/gridpt
                minCost=optVal_initE + c1*E_MAX(1) + c2*E_MAX(2);
                opt_E_mult=[size1_mult,size2_mult];
                optE_SIZE=[E_MAX(1),E_MAX(2)];
            end

            %Store in matrix
            totCost(E1_counter,E2_counter)=(optVal_initE +  c1*E_MAX(1) + c2*E_MAX(2));

            Phi_size{size_iter}=[size(Phi,1),size(Phi,2)];
            
            %end

        end
    end
    
    optRatio(cost_mult_1,cost_mult_2)=optE_SIZE(1)/optE_SIZE(2);
    optCosts{mult_cost_idx}=totCost;
    
    end
end

%Visualize all possible values C1,C2
c1=1:1:size(optRatio,1); c2=1:1:size(optRatio,2);

%Plot optimal ratios vs C1,C2
figure
surf(c2_fact*(c2-1),c1_fact*(c1-1),optRatio(c1,c2));
xlabel('Cost factor c_{2} (1/kWh)'); ylabel('Cost factor c_{1} (1/kWh)'); zlabel('Ratio of sizes, E_{1}^{max}/E_{2}^{max}');
title(sprintf('Variation in optimal sizing with financial costs of storages'));

%OLD: plot optimal cost as a function of storage size
%{
%Visualize all possible policies
max_E1=fliplr(1:1:size(totCost,1)); max_E2=fliplr(1:1:size(totCost,2));

surf((max_E2)*SCALE_SC,(max_E1)*SCALE_BATT,totCost(max_E1,max_E2));
xlabel('Supercapacitor Size (E_2^{max}, kWh)'); ylabel('Battery Size (E_1^{max}, kWh)'); zlabel('Total Cost');
title(sprintf('Optimal Cost as a Function of Storage Size (c1=%d, c2=%d)',c1,c2));
%}