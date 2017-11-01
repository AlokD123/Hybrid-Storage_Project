function [ minCostOut ] = minCostFunc( x,currStatePtr,currStateCostPtr )
%minCostFunc Calculates min cost, all-in-one (version 5)
%   Cost calculation dependent only on x (u set dynamically)
%   Added state 
%   Added pointer passing to avoid repeated recursion
%   ERROR: fminbnd not returning

nextStateCostPtr = RefValue; %Will hold cost of next state
nextStateCostPtr.data = 0; %By default, no cost of next state

nextStatePtr = RefValue;
nextStatePtr.data = currStatePtr.data+1;

f=@(x,u)x+u;

if currStatePtr.data<1
    %In any but last state, find optimum ctrl
    totalCtrlCost=@(u)(minCostFunc5(f(x,u),nextStatePtr,nextStateCostPtr));
    fminbnd(totalCtrlCost,-x,5-x); %<--------- ERROR HERE
else
    %In last state, no control
    uOpt=0;
end

fprintf('State=%d, uOpt=%f',currStatePtr.data,uOpt); %Display state #, optimal ctrl

%Calculate total state cost
currStateCostPtr.data=x^2+uOpt^2+nextStateCostPtr.data;
minCostOut=currStateCostPtr.data;

end
