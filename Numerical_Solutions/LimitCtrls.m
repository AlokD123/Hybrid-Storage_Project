%To limit control (discharge), if would lead to a state out of bounds
function [ D1Opt,D2Opt ] = LimitCtrls( E1,E2,D1Opt,D2Opt,t ) % Input: states, and controls (before saturation), at time t
  global MAX_DISCHARGE; global E_MAX; global E_MIN;
  global expL_State;
  
  %For D1 control...
  if(StateEqn1(E1,D1Opt)>E_MAX(1)) %If next state index would be higher than index of MAX_STATE... (should NOT be possible)
    D1Opt=0;                       %%% NEED TO CHECK!!! %Set control so index of next state is maximum state's index
  elseif(StateEqn1(E1,D1Opt)<E_MIN(1)||D1Opt>MAX_DISCHARGE(1))              %If next state index would be lower than index of MIN_STATE -or- discharge too high
    D1Opt=min( GetCtrl1_CurrNextState(E1,E_MIN(1)), MAX_DISCHARGE(1));      %Set control so index of next state is minimum state's index, or discharge is limited
  end
  %Repeated for D2 control...
  if(StateEqn2(E2,D1Opt,D2Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t))>E_MAX(2))
    D2Opt=GetCtrl2_CurrNextState(E2,E_MAX(2),D1Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t));                  %Set control so index of next state is maximum state's index <------------ %%% NEED TO CHECK that D2>=0 !!!!!
  elseif(StateEqn2(E2,D1Opt,D2Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t))<E_MIN(2)||D2Opt>MAX_DISCHARGE(2))
    D2=min( GetCtrl2_CurrNextState(E2,E_MIN(2),D1Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t)), MAX_DISCHARGE(2)); %Set control so index of next state is minimum state's index, or discharge is limited
    %If still leads to out of bounds, increase nextE2 up to E_MAX(2)-1, and then decrease D1 down to 0
    % (IN THAT ORDER, to minimize cost)
    if( StateEqn2(E2,D1Opt,D2,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t))>E_MAX(2) )
        nextE2=E_MIN(2);
        while(nextE2<=(E_MAX(2)-1) && GetCtrl2_CurrNextState(E2,nextE2,D1Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t))>MAX_DISCHARGE(2))
            nextE2=nextE2+1;
        end
        D1=D1Opt;
        while(D1>=1 && GetCtrl2_CurrNextState(E2,nextE2,D1,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t))>MAX_DISCHARGE(2) && StateEqn1(E1,D1)<E_MAX(1))
            if(  GetCtrl2_CurrNextState(E2,nextE2,D1-1,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t)) >=0) %If discharging remains positive, decrease
                D1=D1-1;
            else
                if(GetCtrl2_CurrNextState(E2,nextE2-1,D1Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t)) < MAX_DISCHARGE(2))
                    nextE2=nextE2-1;
                end
                break;
            end
        end
        D1Opt=D1;
        D2Opt=GetCtrl2_CurrNextState(E2,nextE2,D1Opt,expL_State(E1-E_MIN(1)+1,E2-E_MIN(2)+1,t));
    else
        D2Opt=D2;
    end
  end    
end