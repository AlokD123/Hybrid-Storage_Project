%Storage sizing test with ALP optimization
%---------- Using APPROX LP ONLY!! -----------
%ADDED Regenerative Braking

SCALE_BATT=33.7/250; %0.135 kWh/unit for battery
SCALE_SC=0.16/250; %0.00064 kWh/unit for supercap

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[2,2];
c1=0; c2=200*0;

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

%Approximation sizes
Phi_size={};

global RES_E1; global RES_E2; global RES_L; global RES_U1; 

%for max_E1=2:(1/RES_E1):max_E_SIZE(1)
    %for max_E2=2:(1/RES_E2):max_E_SIZE(2)
        
for max_E1=max_E_SIZE(1):-1:2
    for max_E2=max_E_SIZE(2):-1:2
    %Go through feasible set for E_SIZE
    size1_mult=max_E1; size2_mult=max_E2;
    
    E_MAX=[200*size1_mult;1*size2_mult];
    
    size_iter=size_iter+1; %Next size up
    
    %Define resolutions in simulation
    RES_E1=2/(200*size1_mult);
    RES_E2=2/(size2_mult);
    RES_L=2/(2*size1_mult+3*size2_mult);
    RES_U1=1/size1_mult;
    
    ApproxLP_sol_IHDP_v17_alt; %Get optimal values for this size
    %GetCtrlPolicy_OptQVals_v2; %Get optimal policy matrix
    
    %Store optimal values, for reference;
    optVal_size{size_iter}=ConvCosts; %<--------- DIFFERENCE IS LARGE
    
    %Get optimal value for INITIAL energy state w/ ZERO LOAD
    optVal_initE=ConvCosts(1*(max_E1-E_MIN(1))+1,1*(max_E2-E_MIN(2))+1,(-MIN_LOAD)*RES_L+1);
    
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
    
    Phi_size{size_iter}=[size(Phi,1),size(Phi,2)];
        
    end
end

%Visualize all possible policies
%max_E1=1:(1/RES_E1):max_E_SIZE(1)-1;
%max_E2=1:(1/RES_E2):max_E_SIZE(2)-1;
max_E1=fliplr(1:2:max_E_SIZE(1)-1);
max_E2=fliplr(1:2:max_E_SIZE(2)-1);

%Plot
figure
%Get data differently from matrix depending on whether states are fractional (high resolution) or not
if RES_E1<=1 && RES_E2<=1
    surf((max_E2+1)*SCALE_SC,(max_E1+1)*SCALE_BATT,totCost(max_E1+1,max_E2+1));
elseif RES_E1>1 && RES_E2>1
    surf((max_E2+1)*SCALE_SC,(max_E1+1)*SCALE_BATT,totCost(RES_E1*max_E1,RES_E2*max_E2));
else
    disp('Error!'); %For other cases, break, since not implemented
end

xlabel('Supercapacitor Size (E_2^{max})'); ylabel('Battery Size (E_1^{max})'); zlabel('Total Cost');
title(sprintf('Optimal Cost as a Function of Storage Size (c1=%d, c2=%d)',c1,c2));