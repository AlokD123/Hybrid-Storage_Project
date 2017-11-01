function [ minCost ] = getMinCost( x,k,currIterCostPtr )
%getMinCost Calculates the "cost" metric of an iteration
%   k: iteration #
%   x: "state" value of iteration (updated at each iteration)
%   currIterCostPtr: reference to the cost of the current state
%   minCost: cost of iteration k (output)

%   Note: passing cost of state by reference to avoid repeated recursion
%   (i.e. no need to redo calculating the cost - nextIterCostPtr.data - at
%   the end of the iteration in line 38.)
%   Reason: minimization must be done to obtain values of a "control"
%   variable 'u'. A pointer is the only way to get a return value from the 
%   function 'totalCtrlCost' (assuming the handle doesn't also serve as
%   variable holding the return value?).

%   ERROR: fminbnd not returning???

nextIterCostPtr = RefValue; %nextIterCost holds cost of next iteration
nextIterCostPtr.data = 0;   %Initialize to zero

LAST_ITER=1;    %Recuse for 2 iterations (0 and 1)

if k<LAST_ITER
    %In any but last iteration, find optimum value of u basd on next iteration cost....
    totalCtrlCost=@(u) getMinCost(x+u,k+1,nextIterCostPtr) ; %Get cost of next iteration (k+1). x+u is expression to find value of next state.
    uOpt=fminbnd(totalCtrlCost,-x,5-x);                      %Get optimal value of u, within arbitrary constraints -x<u<-x+5
    %^ ERROR HERE?? (not returning from recursion here)
    y=1; %Does NOT get here
else
    %In last iteration, no "control" needed, so....
    uOpt=0;
end

fprintf('State=%d, uOpt=%f',k,uOpt); %Should display iteration #, optimal u value for the iteration


%Calculate total cost of iteration...
currIterCostPtr.data=x^2+uOpt^2+nextIterCostPtr.data; %"Cost function"
%^ Calculation for cost, in terms of updated state value, optimal u (0 in last iteration), and cost of the following iteration (0 in last iteration)
% NOTE: nextIterCostPtr is local to the iteration, and so should contain
% cost of next iteration after returning from the next iteration at line 26 (not re-initialized)
minCost=currIterCostPtr.data;       %Return value for function (SAME as currIterCostPtr.data)

end