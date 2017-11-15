%warning('off', 'Octave:possible-matlab-short-circuit-operator');
INF=1e6;
clear all;

%Input: initial state
x_init=3; %Must be between MIN_STATE and MAX_STATE

%ASSUMING FINAL COST = 0

%Model setup
E_MIN=[]; %Minimum energy to be stored (lower bound)
E_MAX=[]; %Maximum energy to be stored (upper bound)

MIN_LOAD=0;
MAX_CHARGE=[0;_]; %Maximum charging of the supercapacitor
MAX_DISCHARGE=[]; %Maximum discharging of the 1) battery and 2) supercap
P_PERTURB=normpdf(100,1); %1/(MAX_PERTURB-MIN_PERTURB+1);

j=1; %storage device number (1 or 2)
NUM_STORAGES=2; %number of storages to be used
ALPHA_C=[]; %Efficiency of charging
ALPHA_D=[]; %Efficiency of discharging
BETA=[];    %Storage efficiency


%DP Setup... duplication for each storage
V1(1:(E_MAX(1)-E_MIN(1)+1),1:LAST_ITER) = INF; %COST MATRIX....V(:,k) holds the cost of the kth iteration for each possible state
V2(1:(E_MAX(2)-E_MIN(2)+1),1:LAST_ITER) = INF;

%------------------------STOPPED UPDATING HERE Nov.15,9am ------------------------------------

uOptState(1:(MAX_STATE-MIN_STATE+1),1:LAST_ITER)=0;         %uOptState holds the optimal control for each state, and for all iterations
uOpt(1:LAST_ITER)=0;                            %uOpt holds the optimal control for each iteration, starting from the GIVEN INITIAL STATE
optState(1:LAST_ITER)=INF;                      %holds the value of the optimal state AT a GIVEN iteration

V(:,LAST_ITER)=0; %final cost is 0, for all possible states
%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation

%Define ending iteration
LAST_ITER=4;    %Recurse for 4 iterations (1,2,3,4)



%Set initial conditions
E(1:NUM_STORAGES,1:LAST_ITER)=0;


for t=(LAST_ITER-1):-1:1         %Start at 2nd-last iteration (time, t), and continue backwards
  for state_Index=1:(MAX_STATE-MIN_STATE+1)   %For each state at an iteration...
    x=MIN_STATE+(state_Index-1);      %Map state index to state (TO DO: customize to non-consecutive AND/OR non-integer states)
    CostX(state_Index)=INF;           %CostX will be LOWEST cost-to-go for a state. (Assume infinite cost by default)
    
    for u=(MIN_STATE-MAX_STATE):(MAX_STATE-MIN_STATE)     %For each possible control for that state (at a given iteration)... (TO DO: customize set of possible control)
      CostX_U=0;                      %CostX_U will hold EXPECTED cost of next state, for GIVEN u.
      
      %Find cost-to-go for each possible control, under conditions of perturbations
      %NOTE: this is the perturbation of the next state [w(k+1)]
      for w=MIN_PERTURB:MAX_PERTURB   %Find expected cost-to-go for a given control to be the Expected Cost for over all random perturbances
        if((x+u+w)<=MAX_STATE && (x+u+w)>=MIN_STATE) %If next state is amongst those achievable with a given perturbance....
          nextState_Index=(x+u+w)-MIN_STATE+1;       %Map state to state index, to find cost of next state based on its index
          CostX_U=CostX_U+V(nextState_Index,t+1)*P_PERTURB;   %Add to running cost     (TO DO: customize probability distribution)
        else
          %fprintf('Control not allowed: x=%d,u=%d,x+u=%d\n',x,u,x+u);
        end
      end
      %NOTE: IF NO PERTURBATION.... CostX_U should just hold cost of next state for the given value of u.
      if(CostX_U==0 && (t~=LAST_ITER-1)) %If cannot go to any next state FOR GIVEN CONTROL, update to cost=INF. Exception: in penultimate iteration (ASSUMING final cost=0)
        %ASSUMPTION: all of next possible states have non-zero cost, except for last iteration
        %disp('No achievable next state for this control');
        CostX_U = INF;
      end
      
      
      %Find the Cost-To-Go+Control combination yielding lowest cost for state, amongst those possible for each admissible value of control
      if( (CostX(state_Index)+uOptState(state_Index,t)^2) > (CostX_U+u^2) )  %MOST IMPORTANT
        % (TO DO: customize current control cost calculation)
        CostX(state_Index)=CostX_U;     %For best combo, set cost-to-go to be that cost-to-go for the state
        uOptState(state_Index,t)=u;     %For best combo, set control to give that combo
      end
    end
    
    %Set next state costs in Cost Matrix
    if(CostX(state_Index)==INF)
      disp('State with no achievable next state');
    end    
    V(state_Index,t)=CostX(state_Index) + (x^2+uOptState(state_Index,t)^2);  %Fill in value vector (for a given time t) with: LOWEST costs-to-go for each state + Cost of the state ITSELF (x^2+u^2)
    % (TO DO: customize current state cost calculation)
  end
end


%Get index associated with initial state
initStateInd = x_init-MIN_STATE+1;
%Get minimum total cost given that starting state
NetCost=V(initStateInd,1);
%Get control policy
state_Index=initStateInd; %Starting from initial state...
for(t=1:(LAST_ITER)) %Iterate through cost matrix to find optimal control values
  uOpt(t)=uOptState(state_Index,t);       % Store in uOpt
  optState(t)= MIN_STATE+(state_Index-1); % Also store optimal subsequent states (for reference)
  state_Index = state_Index+uOpt(t);      %(TO DO: ADD EFFECT OF PERTURBATION IN CURRENT STATE)
end



function [ nextE ] = StateEqn( E,u,w,j ) % Input: E(:,t)
  nextE(j)=BETA(j)*E(j)+ALPHA_C(j)*uC(j)-1/ALPHA_D(j)*uD(j)
end