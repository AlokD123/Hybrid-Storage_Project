%Storage sizing test with ALP optimization
%---------- Using APPROX LP ONLY!! -----------
%ADDED Regenerative Braking

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[6,6];
c1=0; c2=0;

INFCOST=1e6;

global E_MAX; global E_MIN;

minCost=Inf; %Default: initial value for J(initial state)+cS cost
vectS_netOptVal=[];%Store CONSTANT state optimal value for each given S
totCost=[]; %Store total cost for each size (convex function)


size_iter=0; %Storage size iteration
PF_opt_mtx=[]; %Transition weights mtx
%optNextE_arr=[]; %Next states array
g_opt_vect=[]; %Stage costs vector
%Store transition weights and stage costs for optimal policies for EACH iteration
PF_opt={}; g_opt={}; Exp_CostToGo={}; optCost_size={}; g_opt_mtx={}; Exp_CostToGo_mtx={};

global RES_E1; RES_E1=1/5;
global RES_E2; RES_E2=1/1000;
global RES_L; RES_L=1;

for max_E1=2:(1/RES_E1):max_E_SIZE(1)
    for max_E2=2:(1/RES_E2):max_E_SIZE(2)
    %Go through feasible set for E_SIZE
    size1_mult=max_E1; size2_mult=max_E2;
    
    E_MAX=[1*size1_mult;1*size2_mult];
    
    size_iter=size_iter+1; %Next size up
    
    ApproxLP_sol_IHDP_v17; %Get optimal values for this size
    %GetCtrlPolicy_OptQVals_v2; %Get optimal policy matrix
    
    %Store optimal values, for reference;
    optVal_size{size_iter}=ConvCosts; %<--------- DIFFERENCE IS LARGE
    
    %Get optimal value for INITIAL energy state w/ ZERO LOAD
    optVal_initE=ConvCosts(RES_E1*(max_E1-E_MIN(1))+1,RES_E2*(max_E2-E_MIN(2))+1,(-MIN_LOAD)*RES_L+1);
    
    vectS_netOptVal=[vectS_netOptVal;optVal_initE]; %Store in vector
    
    %Get optimal storage size till this point
    if (optVal_initE + c1*max_E1 + c2*max_E2) < minCost
        minCost=optVal_initE + c1*max_E1 + c2*max_E2;
        opt_E_SIZE=[max_E1,max_E2];
    end
    
    %Get total cost (convex)
        %Store differently in matrix depending on whether states are fractional (high resolution) or not
        if RES_E1<=1 && RES_E2<=1
            totCost(max_E1,max_E2)=(optVal_initE + c1*max_E1 + c2*max_E2);
        elseif RES_E1>1 && RES_E2>1
            totCost(RES_E1*(max_E1-1),RES_E2*(max_E2-1))=(optVal_initE + c1*max_E1 + c2*max_E2);
        else
            disp('Error!'); %For other cases, break, since not implemented
        end
    
    end
end

%Visualize all possible policies
max_E1=1:(1/RES_E1):max_E_SIZE(1)-1;
max_E2=1:(1/RES_E2):max_E_SIZE(2)-1;

%Plot
figure
%Get data differently from matrix depending on whether states are fractional (high resolution) or not
if RES_E1<=1 && RES_E2<=1
    surf(max_E2,max_E1,totCost(max_E1+1,max_E2+1));
elseif RES_E1>1 && RES_E2>1
    surf(max_E2,max_E1,totCost(RES_E1*max_E1,RES_E2*max_E2));
else
    disp('Error!'); %For other cases, break, since not implemented
end

xlabel('Supercapacitor Size (E_2^{max})'); ylabel('Battery Size (E_1^{max})'); zlabel('Total Cost');
title(sprintf('Optimal Cost as a Function of Storage Size (c1=%d, c2=%d)',c1,c2));