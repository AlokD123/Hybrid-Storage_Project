%warning('off', 'Octave:possible-matlab-short-circuit-operator');

clear all;

%ASSUMING FINAL COST = 0

INF=1e6;
MIN_STATE=-6;
MAX_STATE=6;

MIN_PERTURB=0;
MAX_PERTURB=0;
P_PERTURB=1/(MAX_PERTURB-MIN_PERTURB+1);

%Define ending iteration
LAST_ITER=3;    %Recurse for 3 iterations (1,2,3)

V(1:(MAX_STATE-MIN_STATE+1),1:LAST_ITER) = INF; %V(:,k) holds the cost of the kth iteration for each possible state
uOptState(1:(MAX_STATE-MIN_STATE+1))=0;         %uOptState holds the optimal control for each state, at a SINGLE GIVEN iteration
uOpt(1:LAST_ITER)=0;                            %uOpt holds the optimal control for each iteration
optIterCost(1:LAST_ITER)=INF;                   %holds lowest cost amongst all states AT a GIVEN iteration (i.e. state to choose)
optState(1:LAST_ITER)=INF;                      %holds the value of the optimal state AT a GIVEN iteration

V(:,LAST_ITER)=0;%linspace(1,(MAX_STATE-MIN_STATE+1),13);                               %final cost is least for state -6 //%final cost is 0, for all possible states

%NOTE: at end, uOpt will have best control policy, and optIterCost(1) will contain total cost of DP operation

for t=(LAST_ITER-1):-1:1         %Start at 2nd-last iteration (time, t), and continue backwards
  for state_Index=1:(MAX_STATE-MIN_STATE+1)   %For each state at an iteration...
    x=MIN_STATE+(state_Index-1);      %Map state index to state (TO DO: customize to non-consecutive AND/OR non-integer states)
    CostX(state_Index)=INF;           %CostX will be LOWEST cost-to-go for a state. (Assume infinite cost by default)
    for u=(MIN_STATE-MAX_STATE):(MAX_STATE-MIN_STATE)     %For each possible control for that state (at a given iteration)... (TO DO: customize set of possible control)
      CostX_U=0;                      %CostX_U will hold EXPECTED cost of next state
      for w=MIN_PERTURB:MAX_PERTURB   %Find expected cost-to-go for a given control to be the Expected Cost for over all random perturbances
        if((x+u+w)<=MAX_STATE && (x+u+w)>=MIN_STATE) %If next state is amongst those achievable with a given perturbance....
          nextState_Index=(x+u+w)-MIN_STATE+1;       %Map state to state index, to find cost of next state based on its index
          CostX_U=CostX_U+V(nextState_Index,t+1)*P_PERTURB;   %Add to running cost     (TO DO: customize next-state calculation)
        else
          fprintf('Control not allowed: x=%d,u=%d,x+u=%d\n',x,u,x+u);
        end
      end
      if(CostX_U==0 && (t~=LAST_ITER-1)) %If cannot go to any next state FOR GIVEN CONTROL, update to cost=INF. Exception: in penultimate iteration (ASSUMING final cost=0)
        %ASSUMPTION: all of next possible states have non-zero cost, except for last iteration
        disp('No achievable next state for this control');
        CostX_U = INF;
      end
      if(CostX(state_Index)>CostX_U)    %Find lowest cost-to-go for the state amongst those possible for each of the controls
        CostX(state_Index)=CostX_U;      
        uOptState(state_Index)=u;       %Find best control for the state (providing the lowst cost to go)
      end
    end
    if(CostX(state_Index)==INF)
      disp('State with no achievable next state');
    end    
    V(state_Index,t)=CostX(state_Index) + (x^2+uOptState(state_Index)^2);  %Fill in value vector (for a given time t) with: LOWEST costs-to-go for each state + Cost of the state ITSELF (x^2+u^2)
    % (TO DO: customize current state cost calculation)
    % Possible error: not choosing control to minimize current state cost too (only doing for Cost-to-Go?)
    if(optIterCost(t)>V(state_Index,t)) %Find lowest cost amongst all the states at a given time
      optIterCost(t)= V(state_Index,t);
      uOpt(t)=uOptState(state_Index);   %Optimal control at an iteration amongst all the states is that which has the lowest cost of state.
      optState(t)=MIN_STATE+(state_Index-1);    %Find optimal state at iteration t, based on cost
    end
  end
  %disp('Reached here');
  %disp('Reached here2');
end