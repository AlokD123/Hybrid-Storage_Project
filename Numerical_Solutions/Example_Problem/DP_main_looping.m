%warning('off', 'Octave:possible-matlab-short-circuit-operator');

clear all;

%Input: initial state
x_init=3; %Must be between MIN_STATE and MAX_STATE

%ASSUMING FINAL COST = 0

INF=1e6;
MIN_STATE=1;
MAX_STATE=3;

MIN_PERTURB=0;
MAX_PERTURB=0;
P_PERTURB=1/(MAX_PERTURB-MIN_PERTURB+1);

%Define ending iteration
LAST_ITER=4;    %Recurse for 4 iterations (1,2,3,4)

V(1:(MAX_STATE-MIN_STATE+1),1:LAST_ITER) = INF; %COST MATRIX....V(:,k) holds the cost of the kth iteration for each possible state
uOptState(1:(MAX_STATE-MIN_STATE+1),1:LAST_ITER)=0;         %uOptState holds the optimal control for each state, and for all iterations
uOpt(1:LAST_ITER)=0;                            %uOpt holds the optimal control for each iteration, starting from the GIVEN INITIAL STATE
optState(1:LAST_ITER)=INF;                      %holds the value of the optimal state AT a GIVEN iteration

V(:,LAST_ITER)=0; %final cost is 0, for all possible states
%V(:,LAST_ITER)=linspace(1,(MAX_STATE-MIN_STATE+1),13); %final cost is least for state -6

%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation


for t=(LAST_ITER-1):-1:1         %Start at 2nd-last iteration (time, t), and continue backwards
  for state_Index=1:(MAX_STATE-MIN_STATE+1)   %For each state at an iteration...
    x=MIN_STATE+(state_Index-1);      %Map state index to state (TO DO: customize to non-consecutive AND/OR non-integer states)
    CostX(state_Index)=0;             %CostX will be EXPECTED cost-to-go for a given state.
      
    %Find cost-to-go for each possible value of perturbation (w)
    %NOTE: this is the perturbation of the current time, leading to an expected cost-to-go for the PREV time
    for w=MIN_PERTURB:MAX_PERTURB
      CostX_W=INF;                  %CostX_W will be LOWEST cost of next state, for GIVEN perturbation w. (Assume infinite cost by default)        
      
      for u=(MIN_STATE-MAX_STATE):(MAX_STATE-MIN_STATE) %For each possible control for that state (at a given iteration)...
      %(^TO DO: customize set of possible control)
      
        if((x+u+w)<=MAX_STATE && (x+u+w)>=MIN_STATE) %If next state is amongst those achievable with a given perturbance....
          nextState_Index=(x+u+w)-MIN_STATE+1;       %Map state to state index, to find cost of next state based on its index
          %Find the Cost-To-Go+Control combination yielding lowest cost for state, amongst those possible for each admissible value of control
          if( (CostX_W+uOptState(state_Index,t)^2) > (V(nextState_Index,t+1)+u^2) )  %MOST IMPORTANT.
            %(^TO DO: customize current control cost calculation)
            %Note: Cost-to-go for given u&w (state at next time) is that at time t+1 in matrix V
            
            CostX_W=V(nextState_Index,t+1);     %For best combo, set cost-to-go to be that cost-to-go for the state
            uOptState(state_Index,t)=u;         %For best combo, set control to give that combo
          end
        else
          %fprintf('Control u not allowed, FOR GIVEN w: x=%d,u=%d,x+u=%d\n',x,u,x+u);
        end

      end
      
      %NOTE: IF NO PERTURBATION.... CostX_W should just hold cost of next state for the given value of u.
    
      if(CostX_W==INF) %If cannot go to any next state FOR GIVEN PERTURBATION w...
        %IGNORE possibility of such a perturbation
        %disp('Perturbation w too large. No admissible next state.');
      else
        %ELSE... find expected cost-to-go, to be the Expected Cost for over all random perturbations
        %Find expected value by adding to running cost...
        CostX(state_Index) = CostX(state_Index) + CostX_W*P_PERTURB;
        %(^TO DO: customize probability distribution)
      end
    end
    
    %Set next state costs in Cost Matrix
    if(CostX(state_Index)==INF)
      disp('State with no achievable next state');
    end
    % Cost Function...
    V(state_Index,t)=CostX(state_Index) + (x^2+uOptState(state_Index,t)^2);  %Fill in components of value vector (for a given time t) with: LOWEST cost-to-go for a state + Cost of the state ITSELF (x^2+u^2)
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
  %If control too large/small, limit control...
  if(state_Index+uOpt(t)>(MAX_STATE-MIN_STATE+1)) %If next state index would be higher than index of MAX_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is maximum state's index
  elseif(state_Index+uOpt(t)<(MIN_STATE-MIN_STATE+1)) %If next state index would be lower than index of MIN_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is minimum state's index
  end
  state_Index = state_Index+uOpt(t);
end