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

getMinCost2(5,0);
