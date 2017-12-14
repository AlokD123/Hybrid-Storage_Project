%If limited in control, recalculate resulting optimal subsequent state
function [ nextE1,nextE2 ] = optNextStateLimited( E1,E2,D1Opt,D2Opt,t ) % Input: states, and controls (before saturation), at time t
  global E_MAX;
  global expL_State;
    newE1=StateEqn1(E1,D1Opt);
    if(newE1<0)
        nextE1=0;
    elseif(newE1>E_MAX(1)) %Should not get here
        nextE1=E_MAX(1);
    else
        nextE1=newE1;
    end
    newE2=StateEqn2(E2,D1Opt,D2Opt,expL_State(E1,E2,t));
    if(newE2<0)
        nextE2=0;
    elseif(newE2>E_MAX(2))
        nextE2=E_MAX(2);
    else
        nextE2=newE2;
    end
end