%Calculate subsequent state accounting for no inherent storage loss when
%dropping below lower bound of state
function [ nextE1,nextE2 ] = optNextStateLimited( E1,E2,D1,D2,L ) % Input: state, controls, load
    global E_MIN; global BETA; global ALPHA_D; global ALPHA_C;
    
    newE1=StateEqn1(E1,D1,BETA(1),ALPHA_D(1));
    if(newE1<E_MIN(1))                      %If below lower bound...
        if(StateEqn1(E1,D1,1,ALPHA_D(1))>=E_MIN(1))  %If would not be if with no inherent loss, ASSUME lossless
            nextE1=E_MIN(1);                %In this case, will drop to lower bound
        else
            nextE1=newE1;                   %Else, doesn't matter
        end
    else
        nextE1=newE1;
    end
    
    %Repeat for second storage
    newE2=StateEqn2(E2,D1,D2,L,BETA(2),ALPHA_C(2),ALPHA_D(2));
    if(newE2<E_MIN(2))
        if(StateEqn2(E2,D1,D2,L,1,ALPHA_C(2),ALPHA_D(2))>=E_MIN(2))
            nextE2=E_MIN(2);
        %
        elseif(StateEqn2(E2,D1,D2,L,1,1,1)>=E_MIN(2))&&E2<=2
            nextE2=E_MIN(2);
        %}
        else
            nextE2=newE2;
        end
    else
        nextE2=newE2;
    end
end