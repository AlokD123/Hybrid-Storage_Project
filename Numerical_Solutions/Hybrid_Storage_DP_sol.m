%warning('off', 'Octave:possible-matlab-short-circuit-operator');
clear all;
global INF;
INF=1e6;

global E_MIN; global E_MAX; 
E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
E_MAX=[3;2]; %Maximum energy to be stored (upper bound)

%Input: initial state, horizon
%Initial stored energy (user-defined)
%Must be between MIN_STATE and MAX_STATE
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);
%Recurse for 3 iterations (1,2,3)
LAST_ITER=3;


%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation
%ASSUMING FINAL COST = 0


%Model setup
global MAX_CHARGE; global MAX_DISCHARGE;
MAX_CHARGE=[0;100]; %Maximum charging of the supercapacitor
MAX_DISCHARGE=[2;2]; %Maximum discharging of the 1) battery and 2) supercap

global MIN_LOAD;
MIN_LOAD=0; %Minimum load expected
MAX_LOAD=MAX_DISCHARGE(1)+MAX_DISCHARGE(2);

global ALPHA_C; global ALPHA_D; global BETA; global K;
ALPHA_C=[0.99 0.99]; %Efficiency of charging
ALPHA_D=[0.9;0.95]; %Efficiency of discharging
BETA=[0.99;0.99];    %Storage efficiency
K=2;           %Weighting factor for D1^2 cost
%C1=1;C2=1;     %Cost weighting factors
PERFECT_EFF=0;

%DP Setup... with duplication for each control input
global V; global D1Opt_State; global D2Opt_State; global expCostE;
%COST MATRIX....V(:,k) holds the cost of the kth iteration for each possible state
V(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:LAST_ITER) = INF;       %1 matrix b/c 1 cost function
%uOptState holds the optimal control U for each state, and for all iterations
D1Opt_State(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:LAST_ITER)=0; 
D2Opt_State(1:(E_MAX(1)-E_MIN(1)+1),1:(E_MAX(2)-E_MIN(2)+1),1:(MAX_LOAD-MIN_LOAD+1),1:LAST_ITER)=0;
%final cost is 0, for all possible states and values of "load"
V(:,:,:,LAST_ITER)=0;


%expCostX will be EXPECTED TOTAL for a given state (cost-to-go AND control cost)
%By default, initialize last iteration costs to 0s
expCostE(:,:,LAST_ITER)=V(:,:,1,LAST_ITER);


for t=(LAST_ITER-1):-1:1                %Start at 2nd-last iteration (time, t), and continue backwards
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
        expCostE_L(E_Ind1,E_Ind2,indL)=INF;
        
        %For each possible control for that state (at a given iteration and value of w)...
        %Get CONTROLS and optimal COST of next state (Cost-to-go) for all combos of w and u
        expCostE_L(E_Ind1,E_Ind2,indL)=GetCtrlsUnkNextState( E_Ind1,E_Ind2,indL,t );
                
        %NOTE: IF NO PERTURBATION.... CostX_W should just hold cost of next state for the given value of u.
        if(expCostE_L(E_Ind1,E_Ind2,indL)==INF) %If cannot go to any next state FOR GIVEN PERTURBATION w...
          fprintf('No next state for given L. L=%d, E1=%d, E2=%d\n',L,E_MIN(1)+(E_Ind1-1),E_MIN(2)+(E_Ind2-1));
          %IGNORE possibility of such a perturbation. Perturbation w too large. No admissible next state
        else
          %Else if load admissible, increment count of admissible loads
          numAdmissibleLoads=numAdmissibleLoads+1;
        end
      end
      
      %If the no-load case permits a next state (i.e. not going outside bounds for all controls)...
      if(numAdmissibleLoads~=0)
        %SET PROBABILITY DISTRIBUTION for loads... Uniform  %<------------------------------- **********
        P_PERTURB=1/(numAdmissibleLoads);
        %(^TO DO: customize probability distribution)
        
        %Try to calculate expected cost of the state, now knowing the admissible loads
        for indL=1:(MAX_LOAD-MIN_LOAD+1)
          if(expCostE_L(indL)~=INF) %If CAN go to any next state FOR GIVEN PERTURBATION w...
            %Find expected cost of state, to be the Expected Cost for over all random demands at NEXT TIME STAGE t+1
            %Find expectation by adding to running cost, for each value of load...
            expCostE(E_Ind1,E_Ind2,t) = expCostE(E_Ind1,E_Ind2,t) + V(E_Ind1,E_Ind2,indL,t)*P_PERTURB;
          end
        end
      %Else, if zero-load, zero-control state leads to an expected state out of bounds...
      else
        disp("NO POSSIBLE NEXT STATE FOR CURRENT STATE, for all loads.");
        expCostE(E_Ind1,E_Ind2,t)=INF; %Ignore this state at previous time step when finding the optimal expected next state
      end
      %At end, expCostX contains the expected cost in state (E1,E2)
      
    end
  end
end

%Final costs, depending on load
NetCost=V(E1_INIT-E_MIN(1)+1,E2_INIT-E_MIN(2)+1,:,1);


%uOpt holds the optimal control for each iteration, starting from the GIVEN INITIAL STATE
%D1Opt(1:LAST_ITER)=0;
%D2Opt(1:LAST_ITER)=0;
%optE holds the value of the optimal state (energy stored) AT a GIVEN iteration
%optE1(1:LAST_ITER)=-1;
%optE2(1:LAST_ITER)=-1;


%If zero ("INFINITE") expected cost in given state, next state not possible for given load
%if(CostE1(E_Ind1,E_Ind2)==0)||(CostE2(E_Ind2,E_Ind1)==0 && t~=(LAST_ITER-1)) %ASSUMING NOT POSSIBLE to have zero cost of state for given load EXCEPT FOR 2ND-LAST ITERATION
%  disp('State with no achievable next state for L=%');
%end