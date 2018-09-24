%Storage sizing test with ALP optimization
%---------- Using APPROX LP ONLY!! -----------

%Input: feasible set for E_SIZE, size cost factors
max_E_SIZE=[10,5];
c1=2; c2=1;

global E_MAX; global E_MIN; global MIN_LOAD;

minCost=Inf; %Default: initial value for J(initial state)+cS cost
vectS_netOptVal=[];%Store CONSTANT state optimal value for each given S
totCost=[]; %Store total cost for each size (convex function)


size_iter=0; %Storage size iteration
PF_opt_mtx=[]; %Transition weights mtx
%optNextE_arr=[]; %Next states array
g_opt_vect=[]; %Stage costs vector
%Store transition weights and stage costs for optimal policies for EACH iteration
PF_opt={}; g_opt={}; Exp_CostToGo={}; optCost_size={}; g_opt_mtx={}; Exp_CostToGo_mtx={};


for max_E2=2:max_E_SIZE(2)
    %Go through feasible set for E_SIZE
    max_E1=2*max_E2;            %Choose battery to be at least 2x supercap
    E_MAX=[max_E1;max_E2];
    
    size_iter=size_iter+1; %Next size up
    
    ApproxLP_sol_IHDP_v13; %Get optimal values for this size
    GetCtrlPolicy_OptQVals; %Get optimal policy matrix
    
    %%GET PF and G UNDER OPTIMAL POLICY FOR EACH SIZE
    for l=1:length(E_Ind_VectALL) %For each feasible state...
        %Get TRANSITION WEIGHTS for OPTIMAL POLICY
        VecCtrls=fullPolicyMtx(feasStatesArr(l,1)-E_MIN(1)+1,feasStatesArr(l,2)-E_MIN(2)+1,feasStatesArr(l,3)-MIN_LOAD+1,:);
        %^Get boolean vector of optimal ctrls
        p_opt=find(VecCtrls==1); %Get index of optimal control for state (index p)
        ind=ismember(feasStatesArr(l,:),feasStatesArr_ctrl{p_opt},'rows'); %Get index of state for current control
        PF_opt_mtx=[PF_opt_mtx;full(PF{p_opt}(ind,:))];%Get optimal control transitions weights and store for current state
        %Get NEXT STATES for OPTIMAL POLICY
        %optNextE_arr=[optNextE_arr;full(aug_nextE_Ind_Vect{p_opt}(ind,:))];
        %Get STAGE COSTS for OPTIMAL POLICY
        g_opt_vect=[g_opt_vect;g{p_opt}(ind,:)];
    end
    %Store
    PF_opt{size_iter}=PF_opt_mtx;
    g_opt{size_iter}=g_opt_vect;
    g_opt_mtx{size_iter}=FormatCostVect(g_opt_vect); %<------- Difference is relatively small, in most states
    
    %Reset
    PF_opt_mtx=[]; g_opt_vect=[]; 
    %Get EXPECTED COST-TO-GO
    %Exp_CostToGo{size_iter}=PF_opt{size_iter}*g_opt{size_iter}; <--- NOT expected cost-to-go
    cost_nonneg=cost; cost_nonneg(cost_nonneg<0)=0;
    Exp_CostToGo{size_iter}=PF_opt{size_iter}*cost_nonneg; 
    %Format as matrix
    Exp_CostToGo_mtx{size_iter}=FormatCostVect(Exp_CostToGo{size_iter}); %<------- Difference is bigger
    
    %Store feasible states array, for reference
    feasStatesArr_size{size_iter}=feasStatesArr;
    %Store optimal cost vector, for reference
    optCost_size{size_iter}=cost;
    %Store optimal values, for reference;
    optVal_size{size_iter}=ConvCosts; %<--------- DIFFERENCE IS LARGE
    
    %Get norm of optimal values for CONSTANT energy state
    netOptVal_initE=sum(ConvCosts(size(optVal_size{1},1),size(optVal_size{1},2),1:size(optVal_size{1},3))); %ASSUMING that E_min=0
    
    vectS_netOptVal=[vectS_netOptVal;netOptVal_initE]; %Store in vector
    
    %Get optimal storage size till this point
    if (netOptVal_initE + c1*max_E1 + c2*max_E2) < minCost
        minCost=netOptVal_initE + c1*max_E1 + c2*max_E2;
        opt_E_SIZE=[max_E1,max_E2];
    end
    
    %Get total cost (convex)
    totCost=[totCost;(netOptVal_initE + c1*max_E1 + c2*max_E2)];
    
    %Plot
    plot(max_E2,(netOptVal_initE + c1*max_E1 + c2*max_E2),'O');
    hold on
    
end

title(sprintf('Storage size vs Total Cost (c1=%d, c2=%d)',c1,c2));
xlabel('Size');
ylabel('Total Cost');

%diffVal=optVal_size{4}(1:5,1:3,1:5)-optVal_size{1};