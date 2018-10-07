%Calculate subsequent state accounting for no inherent storage loss when dropping below lower bound of state.
%V2: added regenerative braking
function [ nextE1,nextE2 ] = optNextStateLimited_v2( E1,E2,D1,D2,C2,L ) % Input: state, controls, load
    global E_MIN; global BETA;
    newE1=StateEqn1_wRegenBraking(E1,D1,D2,C2,L,BETA(1));
    if(newE1<E_MIN(1))                                            %If below lower bound...
        if(StateEqn1_wRegenBraking(E1,D1,D2,C2,L,1)>=E_MIN(1))    %If would not be if with no inherent loss, ASSUME lossless
            nextE1=E_MIN(1);                                      %In this case, will drop to lower bound
        else
            nextE1=newE1;                   %Else, doesn't matter
        end
    else
        nextE1=newE1;
    end
    
    %Repeat for second storage
    newE2=StateEqn2_wRegenBraking(E2,D2,C2,BETA(2));
    if(newE2<E_MIN(2))
        if(StateEqn2_wRegenBraking(E2,D2,C2,1)>=E_MIN(2))
            nextE2=E_MIN(2);
        else
            nextE2=newE2;
        end
    else
        nextE2=newE2;
    end
end