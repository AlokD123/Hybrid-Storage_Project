uOptState(1:(MAX_STATE-MIN_STATE+1),1:LAST_ITER)=0;         %uOptState holds the optimal control for each state, and for all iterations
uOpt(1:LAST_ITER)=0;                            %uOpt holds the optimal control for each iteration, starting from the GIVEN INITIAL STATE
optState(1:LAST_ITER)=INF;                      %holds the value of the optimal state AT a GIVEN iteration

V(:,LAST_ITER)=0; %final cost is 0, for all possible states


%If control too large/small, limit control...
  if(state_Index+uOpt(t)>(MAX_STATE-MIN_STATE+1)) %If next state index would be higher than index of MAX_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is maximum state's index
  elseif(state_Index+uOpt(t)<(MIN_STATE-MIN_STATE+1)) %If next state index would be lower than index of MIN_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is minimum state's index
  end