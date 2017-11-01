function [ minCostOut ] = minCostFunc3( x )
%minCostFunc3 Calculates min cost, all-in-one (version 2)
%   Cost calculation dependent only on x (u set dynamically)
%   Issue: repetition of recursion, infinite loop

f=@(x,u)x+u;

totalCtrlCost=@(u)(u^2 + minCostFunc3(f(x,u)));

uOpt=fminbnd(totalCtrlCost,-x,5-x);
fprintf('uOpt=%f',uOpt)

minCostOut=x^2+uOpt^2+ minCostFunc3(f(x,u));

end