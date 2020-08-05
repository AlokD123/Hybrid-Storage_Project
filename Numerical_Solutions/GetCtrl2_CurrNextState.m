%Finds the control value D1 leading to next state E2 from current E2, if BOTH known
%For NO regenerative braking case (uncombined controls)
%Input: E2(t), E2(t+1), D1 (found from D1Opt_State), and L (particular load for which opt control)
function [ D2Opt_State ] = GetCtrl2_CurrNextState( E2,nextE2,D1Opt_State,L )
    global ALPHA_D;global ALPHA_C; global BETA;
    D2=round((nextE2-BETA(2)*E2-ALPHA_C(2)*(D1Opt_State-L))/(ALPHA_C(2)-1/ALPHA_D(2))); %Allow for rounding -0.5->0 up to 0 (to REDUCE discretization issues)
    if D2<0
        disp("Error, D2<0\n");
    end
    D2Opt_State=D2;
end