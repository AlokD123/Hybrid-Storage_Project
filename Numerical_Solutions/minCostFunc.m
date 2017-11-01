function [ minCostOut ] = minCostFunc4( x,k,currStateCostPtr )
%minCostFunc4 Calculates min cost, all-in-one (version 4)
%   Cost calculation dependent only on x (u set dynamically)
%   Added state 
%   Added pointer passing to avoid repeated recursion
%   ERROR: fminbnd not returning

nextStateCostPtr = RefValue; %Will hold cost of next state
nextStateCostPtr.data = 0; %By default, no cost of next state

f=@(x,u)x+u;

if k<3
    %In any but last state, find optimum ctrl
    totalCtrlCost=@(u)(u^2 + minCostFunc4(f(x,u),k+1,nextStateCostPtr));
    uOpt=fminbnd(totalCtrlCost,-x,5-x); %<--------- ERROR HERE
else
    %In last state, no control
    uOpt=0;
end

fprintf('State=%d, uOpt=%f',k,uOpt); %Display state #, optimal ctrl

%Calculate total state cost
currStateCostPtr.data=x^2+uOpt^2+nextStateCostPtr.data;
minCostOut=currStateCostPtr.data;

end