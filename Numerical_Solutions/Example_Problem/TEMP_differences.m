%If control too large/small, limit control...
  if(state_Index+uOpt(t)>(MAX_STATE-MIN_STATE+1)) %If next state index would be higher than index of MAX_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is maximum state's index
  elseif(state_Index+uOpt(t)<(MIN_STATE-MIN_STATE+1)) %If next state index would be lower than index of MIN_STATE...
    uOpt(t) = (MAX_STATE-MIN_STATE+1) - state_Index; %Set control so index of next state is minimum state's index
  end
  
  
  (CostX_W+CtrlCost(D1Opt_State,D2Opt_State,L)) > (V(nextState_Index,t+1)+CtrlCost(D1,D2,L))
  
  
  
  
  
  
%---------------------------------------

%%%Get indiceS associated with initial state
initStateInd = E_INIT-E_MIN+[1 1]';
%%%Get minimum total costS given that starting state
NetCost=[V1(initStateInd(1),1);V2(initStateInd(1),1)];

%Get control policy
state_Ind=initStateInd; %Starting from initial state...
for(t=1:(LAST_ITER)) %Iterate through cost matrices to find optimal control values
  %% Store in uOpt
  D1Opt(t)=D1Opt_State(state_Ind,t);
  D2Opt(t)=D2Opt_State(state_Ind,t);
  %% Also store chain of optimal subsequent states (for reference)
  optE(:,t)= E_MIN+(state_Ind-[1 1]');
  %% Update state index (energy) directly based on the control (discharge)
  %% SUBTRACT because DIScharge
  
  %Implement IF (nextE(j) > BETA(j)*E(j)+ALPHA_C(j)*[uD1+uD2-L]-1/ALPHA_D(j)*uD(j) )--> set uD1&uD2 to LIMIT nextE
  
  state_Ind = state_Ind - [D1Opt(t);D2Opt(t)]; 
  %%%%(TO DO: Update line ↑↑↑ to ADD EFFECT OF PERTURBATION IN CURRENT STATE)
end