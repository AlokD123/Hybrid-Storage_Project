warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear all;
global INF=1e6;

global E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
global E_MAX=[3;1]; %Maximum energy to be stored (upper bound)

%Input: initial state, horizon
%Initial stored energy (user-defined)
%Must be between E_MIN and E_MAX
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);
%Recurse for 3 iterations (1,2,3)
LAST_ITER=3;


%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation
%ASSUMING FINAL COST = 0


%Model setup

global MAX_CHARGE=[0;100]; %Maximum charging of the 1) battery and 2) supercap
global MAX_DISCHARGE=[2;1]; %Maximum discharging of the 1) battery and 2) supercap

MIN_LOAD=0; %Minimum load expected
MAX_LOAD=MAX_DISCHARGE(1)+MAX_DISCHARGE(2); %SHOULD SET MAXIMUM IF DEPENDENT ON STATE?

%j=1; %storage device number (1 or 2)
%NUM_STORAGES=2; %number of storages to be used
global ALPHA_C=[0.99;0.99];%[1;1]; %Efficiency of charging
global ALPHA_D=[0.9;0.95]; %[1;1]; %Efficiency of discharging
global BETA=[0.99;0.99];   %[1;1]; %Storage efficiency
global K=2;           %Weighting factor for D1^2 cost
global C1=1;C2=1;     %Cost weighting factors
PERFECT_EFF=0;

%DP Setup... with duplication for each storage and each control input
global V1; global V2; global D1Opt_State; global D2Opt_State;
%COST MATRIX....V(:,k) holds the cost of the kth iteration for each possible state
V1(1:(E_MAX(1)-E_MIN(1)+1),1:LAST_ITER) = INF; %2 matrices b/c 2 cost functions
V2(1:(E_MAX(2)-E_MIN(2)+1),1:LAST_ITER) = INF;
%uOptState holds the optimal control U for each state, and for all iterations
D1Opt_State(1:(E_MAX(1)-E_MIN(1)+1),1:LAST_ITER)=0; 
D2Opt_State(1:(E_MAX(2)-E_MIN(2)+1),1:LAST_ITER)=0;
%wExpState holds the expected value of the ADMISSIBLE load in each state
LExpState(1:(E_MAX(2)-E_MIN(2)+1),1:LAST_ITER)=0;
%uOpt holds the optimal control for each iteration, starting from the GIVEN INITIAL STATE
D1Opt(1:LAST_ITER)=0;
D2Opt(1:LAST_ITER)=0;
%optE holds the value of the optimal state (energy stored) AT a GIVEN iteration
optE1(1:LAST_ITER)=-1;
optE2(1:LAST_ITER)=-1;
%final cost is 0, for all possible states
V1(:,LAST_ITER)=0;
V2(:,LAST_ITER)=0;


for t=(LAST_ITER-1):-1:1                %Start at 2nd-last iteration (time, t), and continue backwards
  %For each state at an iteration...
  for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
    for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
      %CostX will be EXPECTED cost-to-go for a given state.
      CostE1(E_Ind1)=0;
      CostE2(E_Ind2)=0;
      
      %Reset count of admissible loads
      numAdmissibleLoads=0;
    
      %Find cost-to-go for each possible value of perturbation (w)
      %NOTE: this is the perturbation of the current time, leading to an expected cost-to-go for the PREV time
      for L=MIN_LOAD:MAX_LOAD
        %Map load value to index
        indL=L-MIN_LOAD+1;
        %CostX_W will be LOWEST cost of next state, for GIVEN perturbation w. (Assume infinite cost by default)        
        CostE1_L(indL)=INF;
        CostE2_L(indL)=INF;
        
        %For each possible control for that state (at a given iteration and value of w)...
        %Get optimal cost of next state (and control u leading to that cost) for all combos of w and u
        [CostE1_L(indL),CostE2_L(indL)]=GetCtrlsUnkNextState( E_Ind1,E_Ind2,L,t );
                
        %NOTE: IF NO PERTURBATION.... CostX_W should just hold cost of next state for the given value of u.
        if(CostE1_L(indL)==INF)||(CostE2_L(indL)==INF) %If cannot go to any next state FOR GIVEN PERTURBATION w...
          fprintf('Got here. L=%d, E1=%d, E2=%d\n',L,E_MIN(1)+(E_Ind1-1),E_MIN(2)+(E_Ind2-1));
          %IGNORE possibility of such a perturbation. Perturbation w too large. No admissible next state
        else
          %Else if load admissible, increment count of admissible loads
          numAdmissibleLoads=numAdmissibleLoads+1;
        end
      end
      
      %SET PROBABILITY DISTRIBUTION for loads... Uniform  %<------------------------------- **********
      P_PERTURB=1/(numAdmissibleLoads);
      %(^TO DO: customize probability distribution)
      
      %Try to calculate total expected cost of the state, now knowing the admissible loads
      for indL=1:(MAX_LOAD-MIN_LOAD+1)
        if(CostE1_L(indL)~=INF)&&(CostE2_L(indL)~=INF) %If CAN go to any next state FOR GIVEN PERTURBATION w...
          %Find expected cost-to-go, to be the Expected Cost for over all random perturbations
          %Find expectation by adding to running cost, for each value of load...
          CostE1(E_Ind1) = CostE1(E_Ind1) + CostE1_L(indL)*P_PERTURB;
          CostE2(E_Ind2) = CostE2(E_Ind2) + CostE2_L(indL)*P_PERTURB;
        end
      end
      
      
      %Set next state costs in Cost Matrix
      if(CostE1(E_Ind1)==INF)||(CostE2(E_Ind2)==INF)
        disp('State with no achievable next state');
      end
      %If no single next state (cost) for values of w (i.e. if the expected
      %cost, CostX, falls between possible ones)... Then ROUND to nearest POSSIBLE next cost (state)
      %EXCEPTION: second-last iteration, when all next states have 0 cost...
      if(t~=(LAST_ITER-1))
          %Difference in cost between Expected next cost and Possible costs of next state
          %Default... starting with difference in cost from state 1 as next state
          Diff_Cost1=abs(V1(1,t+1)-CostE1(E_Ind1));
          Diff_Cost2=abs(V2(1,t+1)-CostE2(E_Ind2));
          %Placeholder to hold closest Possible next state cost for given state
          temp_CostE1=V1(1,t+1); 
          temp_CostE2=V2(1,t+1); 
          %Determine actual (possible) control, starting from going to next state index = 1
          D1Opt_State(E_Ind1,t)=E_Ind1-1; %Discharge is subtracted off (must be positive)
          D2Opt_State(E_Ind2,t)=E_Ind2-1;
          %Find closest & lowest cost combination of next states
          for nextE_Ind1=2:(E_MAX(1)-E_MIN(1)+1)
            for nextE_Ind2=2:(E_MAX(2)-E_MIN(2)+1)
              %Find closest-in-cost next states to the expected cost next states
              if abs(V1(nextE_Ind1,t+1)-CostE1(E_Ind1))<Diff_Cost1
                if abs(V2(nextE_Ind2,t+1)-CostE2(E_Ind2))<Diff_Cost2
                  %Find next states minimizing overall cost (WEIGHTED)!!! <------------------ NEED TO CONFIRM!!
                  if((C1*temp_CostE1+C2*temp_CostE2)>(C1*V1(nextE_Ind1,t+1)+C2*V2(nextE_Ind2,t+1)))
                    %Find minimum difference in costs (i.e. to find closest Possible optimal next state)
                    Diff_Cost1=abs(V1(nextE_Ind1,t+1)-CostE1(E_Ind1));
                    Diff_Cost2=abs(V2(nextE_Ind2,t+1)-CostE2(E_Ind2));
                    %Get energies associated with state indices
                    E1=E_MIN(1)+(E_Ind1-1);
                    E2=E_MIN(2)+(E_Ind2-1);
                    nextE1=E_MIN(1)+(nextE_Ind1-1);
                    nextE2=E_MIN(2)+(nextE_Ind2-1);
                    %OBTAIN best possible CONTROL for given state
                    %Pick the higher of the value and zero to ensure POSITIVE
                    D1Opt_State(E_Ind1,t)=max(  round(-ALPHA_D(1)*(nextE1-BETA(1)*E1)) ,  0); %%% NEED TO CHECK!!!
                    D2Opt_State(E_Ind2,t)=max(  round((nextE2-BETA(2)*E2-ALPHA_C(2)*D1Opt_State(E_Ind1,t))/(ALPHA_C(2)-1/ALPHA_D(2))) ,  0);
                    %Get possible next state cost
                    temp_CostE1=V1(nextE_Ind1,t+1);
                    temp_CostE2=V2(nextE_Ind2,t+1);
                  end
                end  
              end
            end
          end
          CostE1(E_Ind1)=temp_CostE1;
          CostE2(E_Ind2)=temp_CostE2;
          %At end, uOptState is a matrix of best controls for each state, and CostX contains lowest possible cost of next state FOR GIVEN STATE
      end
      
      % Cost Function...
      % State Cost = Cost of Next State (lowest cost-to-go&control) + Cost of Control
      %%% TO DO: recheck if L=0 gives cost (I.E. WHETHER would be double counting otherwise??) <-------------------- !!!!!!!!!!!!!!!!!!!!!!!
      V1(E_Ind1,t)=CostE1(E_Ind1) + CtrlCost(D1Opt_State(E_Ind1,t),D2Opt_State(E_Ind2,t),0);
      V2(E_Ind2,t)=CostE2(E_Ind2) + CtrlCost(D1Opt_State(E_Ind1,t),D2Opt_State(E_Ind2,t),0);
      %Fill in components of value vector (for a given time t) with: LOWEST cost of next state for a state + Cost of the control
    end
  end
end

%Get control policy: done in 2 parts...
%Part 1: get state of t=2 (called 'S2') with MINIMUM cost-to-go from GIVEN initial state

%Get index associated with initial state
initE1_Ind = E1_INIT-E_MIN(1)+1;
initE2_Ind = E2_INIT-E_MIN(2)+1;
%Find next state minimizing cost-to-go...
%Start with default values for second state cost and control leading to it
secondCostE1=INF;
secondCostE2=INF;
D1Opt_State(initE1_Ind,1)=0;
D2Opt_State(initE2_Ind,1)=0;
%Variable to hold optimal second state index at end (i.e. index of S2)
optSecondE1_Ind = 0;     
optSecondE2_Ind = 0; 
for secondE1_Ind=1:(E_MAX(1)-E_MIN(1)+1)
  for secondE2_Ind=1:(E_MAX(2)-E_MIN(2)+1)
    %Get energies associated with state indices
    E1=E1_INIT;
    E2=E2_INIT;
    nextE1=E_MIN(1)+(secondE1_Ind-1);
    nextE2=E_MIN(2)+(secondE2_Ind-1);
    
    x=round(-ALPHA_D(1)*(nextE1-BETA(1)*E1)); %%% NEED TO CHECK!!!
    if(PERFECT_EFF==0)
      y=round((nextE2-BETA(2)*E1-ALPHA_C(2)*x)/(ALPHA_C(2)-1/ALPHA_D(2)));
    else
      y=(nextE2-BETA(2)*E1-ALPHA_C(2)*x);
    end
    if( (secondCostE1+CtrlCost(D1Opt_State(initE1_Ind,1),D2Opt_State(initE2_Ind,1),0)) > (V1(secondE1_Ind,2)+CtrlCost(x,y,0)) )  %MOST IMPORTANT. <----------- TO DO: MAKE SAME (nested or not) as "IMPORTANT" above
      if( (secondCostE2+CtrlCost(D1Opt_State(initE1_Ind,1),D2Opt_State(initE2_Ind,1),0)) > (V2(secondE2_Ind,2)+CtrlCost(x,y,0)) )
          %(^TO DO: customize current control cost calculation)
          %Set cost-to-go to be that minimizing cost-to-go for the state
          secondCostE1=V1(secondE1_Ind,2);        
          secondCostE2=V2(secondE2_Ind,2);
          %Set control to lead to optimal second state
          D1Opt_State(initE1_Ind,t)=x;
          D2Opt_State(initE2_Ind,t)=y;
          %Get index of S2
          optSecondE1_Ind=secondE1_Ind;
          optSecondE2_Ind=secondE2_Ind;
      end
    end
  end
end
D1Opt(1)=D1Opt_State(initE1_Ind,1);
D2Opt(1)=D2Opt_State(initE2_Ind,1);
optE1(1)=E1_INIT;
optE2(1)=E2_INIT;
%At end, uOpt(1) contains control leading to optimal second state, and secondCostX contains minimum cost for second state (i.e. state at t=2)


%Part 2: get control policy for remaining times, starting from S2
%Starting from state S2...
E_Ind1=optSecondE1_Ind;
E_Ind2=optSecondE2_Ind;

for(t=2:(LAST_ITER)) %Iterate through cost matrix to find optimal control values
  %Get energies associated with state indices
  E1=E_MIN(1)+(E_Ind1-1);
  E2=E_MIN(2)+(E_Ind2-1);
  % Store in uOpt
  D1Opt(t)=D1Opt_State(E_Ind1,t);
  D2Opt(t)=D2Opt_State(E_Ind2,t);
  % Also store optimal subsequent states (for reference)
  optE1(t)= E1;
  optE2(t)= E2;
  %If control too large/small, limit control...
  %For D1 control...
  if(StateEqn1(E1,D1Opt(t))>E_MAX(1)) %If next state index would be higher than index of MAX_STATE...
    D1Opt(t)=round(-ALPHA_D(1)*(E_MAX(1)-BETA(1)*E1)); %%% NEED TO CHECK!!! %Set control so index of next state is maximum state's index
  elseif(StateEqn1(E1,D1Opt(t))<E_MIN(1))              %If next state index would be lower than index of MIN_STATE...
    D1Opt(t)=round(-ALPHA_D(1)*(E_MAX(1)-BETA(1)*E1)); %%% NEED TO CHECK!!!  %Set control so index of next state is minimum state's index
  end
  %Repeated for D2 control...
  if(StateEqn2(E2,D1Opt(t),D2Opt(t),0)>E_MAX(2))
    D2Opt(t)=round((E_MAX(2)-BETA(2)*E2-ALPHA_C(2)*D1Opt(t))/(ALPHA_C(2)-1/ALPHA_D(2)));
  elseif(StateEqn2(E2,D1Opt(t),D2Opt(t),0)<E_MIN(2))
    D2Opt(t)=round((E_MIN(2)-BETA(2)*E2-ALPHA_C(2)*D1Opt(t))/(ALPHA_C(2)-1/ALPHA_D(2)));
  end
  
  %Update index to next state, after confirming correct control
  %Determine energy of next state...
  nextE1 = StateEqn1(E1,D1Opt(t));
  nextE2 = StateEqn2(E2,D1Opt(t),D2Opt(t),0);
  %Update index to next state...
  E_Ind1 = round(nextE1-E_MIN(1)+1);
  E_Ind2 = round(nextE2-E_MIN(2)+1);
end

%Finally, get minimum total cost given the starting state
%Cost = first state cost + cost to go to second state
NetCost1=CtrlCost(D1Opt(1),D2Opt(1),0)+secondCostE1;
NetCost2=CtrlCost(D1Opt(2),D2Opt(2),0)+secondCostE2;