function [ minCostOut ] = minCostFunc( x,u )
%minCostFunc Calculates min cost

minCostOut=x^2+u^2+minCostFunc(f(x,u), fminbnd(minCostCtrlFunc,-f(x,u),5-f(x,u)) );

end