%IHDP solution (Value Iteration) for Hybrid Storage optimization

%warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear all;

global E_MIN; global E_MAX; 
E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
E_MAX=[10;5]; %Maximum energy to be stored (upper bound)

%Input: initial state, horizon
%Initial stored energy (user-defined)
%Must be between MIN_STATE and MAX_STATE
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);

%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation
%ASSUMING FINAL COST = 0


%Model setup
global MAX_CHARGE; global MAX_DISCHARGE;
MAX_CHARGE=[0;100]; %Maximum charging of the supercapacitor
MAX_DISCHARGE=[10;5]; %Maximum discharging of the 1) battery and 2) supercap

global MIN_LOAD;
MIN_LOAD=0; %Minimum load expected
MAX_LOAD=MAX_DISCHARGE(1)+MAX_DISCHARGE(2);
%SET PROBABILITY DISTRIBUTION for loads... Normal  %<------------------------------- **********
MU_LOAD=floor(0.5*(MAX_LOAD+MIN_LOAD));
%Set stdev so less than 1e-4 probability of outside bounds
SIGMA_LOAD=MAX_LOAD-MIN_LOAD;
while ( (normpdf(MAX_LOAD+1,MU_LOAD,SIGMA_LOAD)>1e-4 || normpdf(MIN_LOAD-1,MU_LOAD,SIGMA_LOAD)>1e-4) && SIGMA_LOAD>1)
    SIGMA_LOAD=SIGMA_LOAD-1;
end
%If probabilities not summing to within 1e-3 of 1, give error
probs=normpdf(linspace(MIN_LOAD,MAX_LOAD,MAX_LOAD-MIN_LOAD+1),MU_LOAD,SIGMA_LOAD);
if sum(probs)<0.999
   disp('Continuous approximation error!!'); 
end
MAX_NUM_ZEROS=3; %Maximum number of zero load counts before end sim

global ALPHA_C; global ALPHA_D; global BETA; global K;
ALPHA_C=[0.99;0.99]; %Efficiency of charging
ALPHA_D=[0.9;0.95]; %Efficiency of discharging
BETA=[0.99;0.99];    %Storage efficiency
K=2;           %Weighting factor for D1^2 cost
PERFECT_EFF=0;
%Recurse for <=MAX_ITER iterations, even if not reached stopping condition for VI
MAX_ITER=10;

%DP Setup... with duplication for each control input
global V; global D1Opt_State; global D2Opt_State; global expCostE;
%COST MATRIX....V(:,k) holds the cost of the kth iteration for each possible state
V(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:MAX_ITER) = Inf;       %1 matrix b/c 1 cost function
%uOptState holds the optimal control U for each state, and for all iterations
D1Opt_State(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:MAX_ITER)=0; 
D2Opt_State(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:MAX_ITER)=0;
%final cost is 0, for all possible states and values of "load"
V(:,:,:,MAX_ITER)=0;

%optNextE will hold optimal NEXT state at state E with load L (at iteration t)... FOR REFERENCE
global optNextE1; global optNextE2;
optNextE1(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:MAX_ITER)=Inf;
optNextE2(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:MAX_ITER)=Inf;

%expCostX will be EXPECTED TOTAL for a given state (cost-to-go AND control cost)
%By default, initialize last iteration costs to 0s
expCostE(:,:,MAX_ITER)=V(:,:,1,MAX_ITER);


%IHDP w/ VALUE ITERATION
%Stopping condition...(Cost of current state)-(Cost of next state) <= VI_ERR
VI_ERR=1;
%Discounted infinite horizon problem
global DISCOUNT; %Discount factor
global BOOL_VI_CONV; %Array of booleans indicating convergence for that state (combination of components)
DISCOUNT=0.99;
BOOL_VI_CONV(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1))=0;
BOOL_VI_CONV_PREV=BOOL_VI_CONV; %holder variable for previous array, to check if sparsity decreasing
%Store difference in cost between iterations in a matrix
diffV(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1))=0;

%STEP 1: Obtain optimal policy for infinite horizon case (OFFLINE)
t=MAX_ITER-1; %Start at 2nd-last iteration (time, t)
while (  t>0 && ~all(all(all(BOOL_VI_CONV(:,:,:)==1))) )               %Continue backwards until VI converges or reach t=0, whichever first
  BOOL_VI_CONV_PREV=BOOL_VI_CONV; %Store previous iteration of VI convergence array
  
  %For each state at an iteration...
  for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
    for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
      %Since expected cost found from adding to running cost, set to 0s initially
      expCostE(E_Ind1,E_Ind2,t)=0;
      %Reset count of admissible loads
      numAdmissibleLoads=0;
    
      %Find cost-to-go for each possible value of perturbation (w)
      %NOTE: this is the perturbation of the current time, leading to an expected cost-to-go for the PREV time
      for indL=1:(MAX_LOAD-MIN_LOAD+1)
        %Map index to value of load
        L=indL+MIN_LOAD-1;
        %CostX_W will be LOWEST cost of next state, for GIVEN perturbation w. (Assume infinite cost by default)        
        expCostE_L(E_Ind1,E_Ind2,indL)=Inf;
        
        %For each possible control for that state (at a given iteration and value of w)...
        %Get CONTROLS and optimal COST of next state (Cost-to-go) for all combos of w and u
        expCostE_L(E_Ind1,E_Ind2,indL)=GetCtrlsUnkNextState( E_Ind1,E_Ind2,indL,t );

        %NOTE: IF NO PERTURBATION.... CostX_W should just hold cost of next state for the given value of u.
        if(expCostE_L(E_Ind1,E_Ind2,indL)==Inf) %If cannot go to any next state FOR GIVEN PERTURBATION w...
          %fprintf('No next state for given L. L=%d, E1=%d, E2=%d\n',L,E_MIN(1)+(E_Ind1-1),E_MIN(2)+(E_Ind2-1));
          %IGNORE possibility of such a perturbation. Perturbation w too large. No admissible next state
        else
          %Else if load admissible, increment count of admissible loads
          numAdmissibleLoads=numAdmissibleLoads+1;
          
          %VALUE ITERATION TEST
          %Check for change in cost for state within VI_ERR
          if((V(E_Ind1,E_Ind2,indL,t)-V(E_Ind1,E_Ind2,indL,t+1))<=VI_ERR) %If change is within error...
            BOOL_VI_CONV(E_Ind1,E_Ind2,indL)=1;
          end
          %Store difference in cost in matrix:
          diffV(E_Ind1,E_Ind2,indL)=V(E_Ind1,E_Ind2,indL,t)-V(E_Ind1,E_Ind2,indL,t+1);
        end
        
        %VI TEST
        %If current or next state cost is INF (if next is not last), ignore change in this state for the test
        if(V(E_Ind1,E_Ind2,indL,t)==Inf) %&& (t+1)~=MAX_ITER)
            BOOL_VI_CONV(E_Ind1,E_Ind2,indL)=1;
        elseif (V(E_Ind1,E_Ind2,indL,t+1)==Inf)
            BOOL_VI_CONV(E_Ind1,E_Ind2,indL)=1;
        end
      end
      
      %If the no-load case permits a next state (i.e. not going outside bounds for all controls)...
      if(numAdmissibleLoads~=0)
        P_PERTURB=1/(numAdmissibleLoads); %Set load distribution to be uniform, by default
        
        %Try to calculate expected cost of the state, now knowing the admissible loads
        for indL=1:(MAX_LOAD-MIN_LOAD+1)
          if(expCostE_L(E_Ind1,E_Ind2,indL)~=Inf) %If CAN go to any next state FOR GIVEN PERTURBATION w...
            %Find expected cost of state, to be the Expected Cost for over all random demands at NEXT TIME STAGE t+1
            %1) Determine probability of given load value
                %Map index to value of load
                L=indL+MIN_LOAD-1;
            %P_PERTURB=normpdf(L,MU_LOAD,SIGMA_LOAD);
            %2) Find expectation by adding to running cost, for each value of load...
            expCostE(E_Ind1,E_Ind2,t) = expCostE(E_Ind1,E_Ind2,t) + V(E_Ind1,E_Ind2,indL,t)*P_PERTURB;
          end
        end
      %Else, if zero-load, zero-control state leads to an expected state out of bounds...
      else
        disp("NO POSSIBLE NEXT STATE FOR CURRENT STATE, for all loads.");
        expCostE(E_Ind1,E_Ind2,t)=Inf; %Ignore this state at previous time step when finding the optimal expected next state
      end
      %At end, expCostX contains the expected cost in state (E1,E2)
      
    end
  end
  
  %VisualizeBool_VI_CONV;
  %VisualizeOptNextState;

  %Visualize the convergence by decrease in matrix norm of 3D difference matrix
  norm_array=arrayfun(@(idx) norm(diffV(:,:,idx)), 1:size(diffV,3));
  norm_diffV(MAX_ITER-t)=norm(norm_array);
  
  %If closer to convergence...
  if( ~all(all(all( BOOL_VI_CONV(:,:,:)>BOOL_VI_CONV_PREV(:,:,:) ))) )
      fprintf("Closer to convergence @t=%d\n",t);
  end
  t=t-1; %Continue bkwds in recursion;
end
%Replace infinite costs with -1%
V(V==inf)=-1;
%Final costs, depending on load
NetCost=V(E1_INIT-E_MIN(1)+1,E2_INIT-E_MIN(2)+1,:,1);

%GET INFINITE HORIZON POLICY
D1Opt_Inf=D1Opt_State(:,:,:,t+1);
D2Opt_Inf=D2Opt_State(:,:,:,t+1);

%Get offset time index (time index when value iteration converged)
IND_T_OFFS=t;
%Restart at last time
t=t+1;


%%STEP 2: evaluate infinite horizon policy for a random sequence of loads (ONLINE)
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads
while t<=(MAX_ITER-1)
    %Set 1-indexed time, with no offset due to VI-stopping
    t_ind_VI=t-(IND_T_OFFS);
    %Set state index
    indE1=optE1(t_ind_VI)-E_MIN(1)+1;
    indE2=optE2(t_ind_VI)-E_MIN(2)+1;
    %Create random demand from IID Uniform probability sequence
    MAX_LOAD_STATE=optE1(t_ind_VI)+optE2(t_ind_VI)-1; %Maximum possible load limited to total energy stored in that state
    if(MAX_LOAD_STATE==Inf)
        MAX_LOAD_STATE=MAX_LOAD;
    end
    L=randi(MAX_LOAD_STATE-MIN_LOAD+1,1,1)+MIN_LOAD-1;
    %Short form for optimal controls...
    D1=D1Opt_Inf(indE1,indE2,L-MIN_LOAD+1);
    D2=D2Opt_Inf(indE1,indE2,L-MIN_LOAD+1);
    %Calculate the state these values of u and w will lead to, even if impossible...
        [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
        
    %If leads to next state guaranteed out of bounds, decrease load
    while(nextE1>E_MAX(1)||nextE1<E_MIN(1)||nextE2>E_MAX(2)||nextE2<E_MIN(2)||(D1+D2-L)<0)
        L=L-1;%Decrement load
        %Increment out of bounds count
        countOOB=countOOB+1;
        if(L==-1)
            L=0;
            break;
        end
        
        %Calculate next state for new load value
        D1=D1Opt_Inf(indE1,indE2,L-MIN_LOAD+1);
        D2=D2Opt_Inf(indE1,indE2,L-MIN_LOAD+1);
        [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    indL=L-MIN_LOAD+1;
    Load(t_ind_VI)=L;          %Hold value of load (for reference)
    %Get optimal costs for this sequence of loads
    optV(t_ind_VI)=V(indE1,indE2,indL,t);
    %Get costs of other possible states for this state and load??
    %(For comparison)

    %Round next state to nearest int
    nextE1=round(nextE1);
    nextE2=round(nextE2);
    %TO DO: APPLY STATE INTERPOLATION for better accuracy
    
    %Get specific control sequence for given load sequence
    D1Opt(t_ind_VI)=D1;
    D2Opt(t_ind_VI)=D2;
    %Get optimal next state, starting at E_INIT
    optE1(t_ind_VI+1)=nextE1;
    optE2(t_ind_VI+1)=nextE2;
    %Note: cross-check with optNextStateLtd?
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        optE1(t_ind_VI+1)=optE1(t_ind_VI); optE2(t_ind_VI+1)=optE2(t_ind_VI);
        optV(t_ind_VI+1)=V(optE1(t_ind_VI+1)-E_MIN(1)+1,optE2(t_ind_VI+1)-E_MIN(2)+1,indL,t);
        Load(t_ind_VI+1)=0;
        t=MAX_ITER;
    end
    t=t+1;
end