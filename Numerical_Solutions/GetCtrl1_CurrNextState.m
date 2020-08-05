%Finds the control value D1 leading to next state E1 from current E1, if BOTH known
%For NO regenerative braking case (uncombined controls)
%Input: E1(t), E1(t+1)
function [ D1Opt_State ] = GetCtrl1_CurrNextState( E1,nextE1 )
    global ALPHA_D; global BETA;
    D1=round(-ALPHA_D(1)*(nextE1-BETA(1)*E1)); %Allow for rounding -0.5->0 up to 0 (to REDUCE discretization issues)
    if D1<0
        disp("Error, D1<0\n");
    end
    D1Opt_State=D1;
end