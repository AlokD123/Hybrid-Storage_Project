clear all;

INF=1e6;
MIN_STATE=-10;
MAX_STATE=10;

MIN_PERTURB=-2;
MAX_PERTURB=2;
P_PERTURB=1/(MAX_PERTURB-MIN_PERTURB+1);

j=1; %storage device number
ALPHA=[];
BETA=[];

global uOpt;
global LAST_ITER;

%Define ending iteration
LAST_ITER=2;    %Recurse for 2 iterations (1 and 2)

V(MIN_STATE:MAX_STATE,1:LAST_ITER) = INF;     %V(:,k) holds the cost of the kth iteration for each possible state
uOptState(1:MAX_STATE-MIN_STATE+1)=0;         %uOptState holds the optimal control for each state, at a SINGLE GIVEN iteration
uOpt(1:LAST_ITER)=0;                          %uOpt holds the optimal control for each iteration
optIterCost(1:LAST_ITER)=INF;                 %holds lowest cost amongst all states AT a GIVEN iteration (i.e. state to choose)

V(:,LAST_ITER)=0;                             %final cost is 0, for all possible states

for t=(LAST_ITER-1):1         %Start at 2nd-last iteration (time, t), and continue backwards
  for x=MIN_STATE:MAX_STATE   %For each state at an iteration...
    CostX(x)=INF;             %CostX will LOWEST cost-to-go for a state. (Assume infinite cost by default)
    for u=-x:(-x+5)           %For each possible control for that state (at a given iteration)... (TO DO: customize set of possible control)
      CostX_U=0;              %CostX_U will hold EXPECTED cost of next state
      for w=MIN_PERTURB:MAX_PERTURB     %Find expected cost-to-go for a given control to be the Expected Cost for over all random perturbances
        if((x+u+w)<MAX_STATE && (x+u+w)>MIN_STATE) %If next state is amongst those achievable with a given perturbance....
          CostX_U=CostX_U+V(t+1,x+u+w)*P_PERTURB   %Add to running cost     (TO DO: customize next-state calculation)
        end
      end
      if(CostX_U=0) disp('State with no achievable next state') return end; %If cannot go to any next state, break script. (TO DO: update to cost=INF)
      if(CostX(x)>CostX_U)    %Find lowest cost-to-go for the state amongst those possible for each of the controls
        CostX(x)=CostX_U;
        uOptState(x)=u;       %Find best control for the state (providing the lowst cost to go)
      end
    end
    V(t,x)=CostX(x) + (x^2+uOptState(x)^2);  %Fill in value vector (for a given time t) with: LOWEST costs-to-go for each state + Cost of the state ITSELF (x^2+u^2)
    % (TO DO: customize current state cost calculation)
    if(optIterCost(t)>V(t,x)) %Find lowest cost amongst all the states at a given time
      optIterCost(t)= V(t,x);
      uOpt(t)=uOptState(x);   %Optimal control at an iteration amongst all the states is that which has the lowest cost of state.
  end
end

function [ nextX ] = StateEqn( x,u,w,1 )
  nextX(1)=BETA(1)*x(1)+ALPHA_C(1)*uC(1)-1/ALPHA_D(1)*uD(1)
end

%NOTE: at end, uOpt will have best control policy, and optIterCost(1) will contain total cost of DP operation