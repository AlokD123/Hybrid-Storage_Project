clear all;

global iterCost;
global uOpt;
global LAST_ITER;

%Define ending iteration
LAST_ITER=1;    %Recurse for 2 iterations (0 and 1)

 for i=0:LAST_ITER
     iterCost=[iterCost;0]; %iterCost(k) holds cost of the kth iteration
     uOpt=[uOpt;0];         %uOpt(k) holds optimal control of the kth iteration
 end

 
 function [ minCost ] = getMinCost2( x,k )
%getMinCost2 Calculates the "cost" metric of an iteration
%   k: iteration #
%   x: "state" value of iteration (updated at each iteration)
%   iterCost: cost of the current iteration
%   minCost: cost of iteration k (output)

global iterCost;
global uOpt;

global LAST_ITER;

if k<LAST_ITER
    %In any but last iteration, find optimum value of u basd on next iteration cost....
    totalCtrlCost=@(u) getMinCost2(x+u,k+1) ; %Get cost of next iteration (k+1). x+u is expression to find value of next state.
    uOpt(k+1)=fminbnd(totalCtrlCost,-x,5-x);       %Get optimal value of u, within arbitrary constraints -x<u<-x+5
    %Calculate total cost of iteration...
    iterCost(k+1)=x^2+uOpt(k+1)^2+iterCost(k+2); %"Cost function"
else
    %In last iteration, no "control" needed, so....
    uOpt(k+1)=2;
    %Calculate total cost of iteration...
    iterCost(k+1)=x^2; %"Cost function"
end

fprintf('State=%d, uOpt=%f, iterCost=%f\n',k,uOpt(k+1),iterCost(k+1)); %Should display iteration #, optimal u value, and cost

%NOTE: Calculation for cost is in terms of updated state value, optimal u (0 in
%last iteration), and cost of the following iteration (0 in last iteration)
minCost=iterCost(k+1);       %Return value for function (SAME as iterCost)

end