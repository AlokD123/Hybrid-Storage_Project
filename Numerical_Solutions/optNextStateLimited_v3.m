%Calculate subsequent state accounting for no inherent storage loss when dropping below lower bound of state.
%V3: single control, with regenerative braking
function [ nextE1,nextE2 ] = optNextStateLimited_v3( E1,E2,U1,L ) % Input: state, control, load
    global E_MIN; global BETA; global ALPHA_C; global ALPHA_D;
    newE1=StateEqn1_wRegenBraking_v2(E1,U1,BETA(1));
    if(newE1<E_MIN(1))                                      %If below lower bound...
        if(StateEqn1_wRegenBraking_v2(E1,U1,1)>=E_MIN(1))    %If would not be if with no inherent loss, ASSUME lossless
            nextE1=E_MIN(1);                                %In this case, will drop to lower bound
        else
            nextE1=newE1;                   %Else, doesn't matter
        end
    else
        nextE1=newE1;
    end
    
    %Repeat for second storage
    newE2=StateEqn2_wRegenBraking_v2(E2,U1,L,BETA(2),ALPHA_C(2),ALPHA_D(2));
    if(newE2<E_MIN(2))
        if(StateEqn2_wRegenBraking_v2(E2,U1,L,1,ALPHA_C(2),ALPHA_D(2))>=E_MIN(2))
            nextE2=E_MIN(2);
        %{
        elseif(StateEqn2(E2,U1,L,1,1,1)>=E_MIN(2))&&E2<=2
            nextE2=E_MIN(2);
        %}
        else
            nextE2=newE2;
        end
    else
        nextE2=newE2;
    end
end