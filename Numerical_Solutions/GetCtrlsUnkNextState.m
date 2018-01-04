function [ expCostE_L ] = GetCtrlsUnkNextState( E_Ind1,E_Ind2,indL,t ) % Input: E1(t),E2(t), and L (particular load for which opt control)
  global MAX_DISCHARGE;global MAX_CHARGE; global E_MAX; global E_MIN; global INF; global MIN_LOAD;
  global V; global D1Opt_State; global D2Opt_State; global expCostE;
  
    %Map state index to state    %%(TO DO: customize to non-consecutive AND/OR non-integer states)
    E1=E_MIN(1)+(E_Ind1-1);
    E2=E_MIN(2)+(E_Ind2-1); 
    %Map index to value of load
    L=indL+MIN_LOAD-1; 
  
    %tempCostX_W will contain LOWEST cost of next state, for GIVEN perturbation w. (Assume infinite cost by default)        
    tempCostE_L=INF;

    %For each possible control for that state (at a given iteration and value of w)...
    for D1=0:MAX_DISCHARGE(1)
      for D2=0:MAX_DISCHARGE(2)
        %Calculate the state these values of u and w will lead to, even if
        %impossible...
        [nextE1,nextE2]=optNextStateLimited(E1,E2,D1,D2,L);
        %If next state is amongst those achievable with a given perturbance....
        if(nextE1<=E_MAX(1) && nextE1>=E_MIN(1))
          if(nextE2<=E_MAX(2) && nextE2>=E_MIN(2))

            %IF meeting following conditions: (C_MAX and C_MIN)
            %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
            if(~((D1+D2-L)<0||(D1+D2-L)>MAX_CHARGE(2)))
              %Map state to state index, to find cost of next state based on its index
              nextE_Ind1=round(nextE1-E_MIN(1)+1);
              nextE_Ind2=round(nextE2-E_MIN(2)+1); 

              %Find the Cost-To-Go+Control combination yielding lowest cost for state, amongst those possible for each admissible value of control
              %MOST IMPORTANT: minimization
              if( (tempCostE_L+CtrlCost(D1Opt_State(E_Ind1,E_Ind2,indL,t),D2Opt_State(E_Ind1,E_Ind2,indL,t),L)) > (expCostE(nextE_Ind1,nextE_Ind2,t+1)+CtrlCost(D1,D2,L)) )
                %Note: Cost-to-go for given u&w (state at next time) is that at time t+1 in matrix V

                %For best combo, set control to give that combo
                D1Opt_State(E_Ind1,E_Ind2,indL,t)=D1;  
                D2Opt_State(E_Ind1,E_Ind2,indL,t)=D2;  
                
                %For best combo, set cost-to-go to be that cost-to-go for the state
                tempCostE_L=expCostE(nextE_Ind1,nextE_Ind2,t+1);
                
                %Set the lowest TOTAL cost of the current state, for given controls
                V(E_Ind1,E_Ind2,indL,t)=tempCostE_L+CtrlCost(D1,D2,L);           %<------ TOTAL COST
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
    expCostE_L=tempCostE_L; %Return optimal cost
    % At end...
    % - expCostX_W contains optimal next state cost, for given load and state
    % - V contains the TOTAL cost of the current state, for given load
    % - uOpt_State contains optimal control policy
end