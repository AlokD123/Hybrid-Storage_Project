%Storage sizing test with ALP optimization
%---------- Using APPROX LP ONLY!! -----------

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[8,8];
c1=0.025; c2=0.25;

global E_MAX; %global E_MIN; global MIN_LOAD;

minCost=Inf; %Default: initial value for J(initial state)+cS cost
vectS_netOptVal=[];%Store CONSTANT state optimal value for each given S
totCost=[]; %Store total cost for each size (convex function)


size_iter=0; %Storage size iteration
PF_opt_mtx=[]; %Transition weights mtx
%optNextE_arr=[]; %Next states array
g_opt_vect=[]; %Stage costs vector
%Store transition weights and stage costs for optimal policies for EACH iteration
PF_opt={}; g_opt={}; Exp_CostToGo={}; optCost_size={}; g_opt_mtx={}; Exp_CostToGo_mtx={};


for max_E1=4:max_E_SIZE(1)
    for max_E2=4:max_E_SIZE(2)
       %Go through feasible set for E_SIZE
        E_MAX=[max_E1;max_E2];

        size_iter=size_iter+1; %Next size up

        ApproxLP_sol_IHDP_v13; %Get optimal values for this size
        GetCtrlPolicy_OptQVals; %Get optimal policy matrix

        %Get optimal values
        optVal_size{size_iter}=ConvCosts;

        %Get norm of optimal values for CONSTANT energy state
        netOptVal_initE=sum(ConvCosts(size(optVal_size{1},1),size(optVal_size{1},2),1:size(optVal_size{1},3))); %ASSUMING that E_min=0

        vectS_netOptVal=[vectS_netOptVal;netOptVal_initE]; %Store in vector

        %Get optimal storage size till this point
        if (netOptVal_initE + c1*max_E1 + c2*max_E2) < minCost
            minCost=netOptVal_initE + c1*max_E1 + c2*max_E2;
            opt_E_SIZE=[max_E1,max_E2];
        end

        %Store total cost (convex) in table for each storage size
        totCost=[totCost; max_E1, max_E2,(netOptVal_initE + c1*max_E1 + c2*max_E2)];

        %Plot
        %plot(max_E2,(netOptVal_initE + c1*max_E1 + c2*max_E2),'O');
        %hold on 
    end
end

%title(sprintf('Storage size vs Total Cost (c1=%d, c2=%d)',c1,c2));
%xlabel('Size');
%ylabel('Total Cost');

%diffVal=optVal_size{4}(1:5,1:3,1:5)-optVal_size{1};