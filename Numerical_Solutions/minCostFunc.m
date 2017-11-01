function [ minCostOut ] = minCostFunc2( x,u )
%minCostFunc2 Calculates min cost, all-in-one

totalCtrlCost=@(u)(u^2 + minCostFunc2(  f(x,u), fminbnd(nextStateCost,-f(x,u),5-f(x,u)) ));
nextStateCost=@(v) minCostFunc2(f(x,u),v);

u0=fminbnd(totalCtrlCost,-x,5-x);

minCostOut=x^2+u0^2+ minCostFunc2(  f(x,u), fminbnd(nextStateCost,-f(x,u),5-f(x,u)) );

end