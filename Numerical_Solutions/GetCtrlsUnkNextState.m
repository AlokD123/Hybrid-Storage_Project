function [ CostE1_L,CostE2_L ] = GetCtrlsUnkNextState( E_Ind1,E_Ind2,L,t ) % Input: E1(t),E2(t), and L (particular load for which opt control)
  global MAX_DISCHARGE;global MAX_CHARGE; global E_MAX; global E_MIN; global INF;
  global V1; global V2; global D1Opt_State; global D2Opt_State;
      %Map state index to state    %%(TO DO: customize to non-consecutive AND/OR non-integer states)
      E1=E_MIN(1)+(E_Ind1-1);
      E2=E_MIN(2)+(E_Ind2-1);  
  
    %CostX_W will be LOWEST cost of next state, for GIVEN perturbation w. (Assume infinite cost by default)        
    CostE1_L=INF;
    CostE2_L=INF;

    %For each possible control for that state (at a given iteration and value of w)...
    for D1=0:MAX_DISCHARGE(1)
      for D2=0:MAX_DISCHARGE(2)
        %If next state is amongst those achievable with a given perturbance....
        if(StateEqn1(E1,D1)<=E_MAX(1) && StateEqn1(E1,D1)>=E_MIN(1))
          if(StateEqn2(E2,D1,D2,L)<=E_MAX(2) && StateEqn2(E2,D1,D2,L)>=E_MIN(2))

            %IF meeting following conditions: (C_MAX and C_MIN)
            %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
            if(~((D1+D2-L)<0||(D1+D2-L)>MAX_CHARGE(2)))
              %Map state to state index, to find cost of next state based on its index
              nextE_Ind1=round(StateEqn1(E1,D1)-E_MIN(1)+1);
              nextE_Ind2=round(StateEqn2(E2,D1,D2,L)-E_MIN(2)+1); 

              %Find the Cost-To-Go+Control combination yielding lowest cost for state, amongst those possible for each admissible value of control
              %MOST IMPORTANT: minimization
              if( (CostE1_L+CtrlCost(D1Opt_State(E_Ind1,E_Ind2,t),D2Opt_State(E_Ind2,E_Ind1,t),L)) > (V1(nextE_Ind1,nextE_Ind2,t+1)+CtrlCost(D1,D2,L)) )
                if( (CostE2_L+CtrlCost(D1Opt_State(E_Ind1,E_Ind2,t),D2Opt_State(E_Ind2,E_Ind1,t),L)) > (V2(nextE_Ind2,nextE_Ind1,t+1)+CtrlCost(D1,D2,L)) )
                  %Note: Cost-to-go for given u&w (state at next time) is that at time t+1 in matrix V

                  %For best combo, set cost-to-go to be that cost-to-go for the state
                  CostE1_L=V1(nextE_Ind1,nextE_Ind2,t+1); 
                  CostE2_L=V2(nextE_Ind2,nextE_Ind1,t+1); 
                  %For best combo, set control to give that combo
                  D1Opt_State(E_Ind1,E_Ind2,t)=D1;  
                  D2Opt_State(E_Ind2,E_Ind1,t)=D2;  
                  %^^ value of uOptState obtained here is UNUSED after value of CostX_W is found!!!
                end
              end
            end
          else
            %fprintf('nextE_2=%d not allowed\n',round(StateEqn2(E2,D1,D2,L))); %Control combo D1,D2 not allowed, FOR GIVEN L
          end
        else
          %fprintf('nextE_1=%d not allowed\n',round(StateEqn1(E1,D1))); %Control D1 not allowed, FOR GIVEN E2
        end
      end
    end
  
end