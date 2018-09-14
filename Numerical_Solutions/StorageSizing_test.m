%Storage sizing test with ALP optimization

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[10,5];
c1=1; c2=1;

global E_MAX;

minCost=Inf; %Default: initial value for J(initial state)+cS cost
vectS_netOptVal=[];%Store initial state optimal value for each given S


for max_E2=3:max_E_SIZE(2)
    %Go through feasible set for E_SIZE
    max_E1=2*max_E2;            %Choose battery to be at least 2x supercap
    E_MAX=[max_E1;max_E2];
    
    ApproxLP_sol_IHDP_v13; %Get optimal values for this size
    
    %Get norm of optimal values for initial (maximum) energy state
    netOptVal_initE=sum(ConvCosts(max_E1+1,max_E2+1,:)); %ASSUMING that E_min=0
    
    vectS_netOptVal=[vectS_netOptVal;netOptVal_initE]; %Store in vector
    
    if (netOptVal_initE + c1*max_E1 + c2*max_E2) < minCost
        minCost=netOptVal_initE;
        opt_E_SIZE=[max_E1,max_E2];
    end
    
    i=i+1;
end