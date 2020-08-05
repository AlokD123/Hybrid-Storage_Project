%Expression that is constant for a given state aggregation. Hand-engineered; can be modified here

function [fitExpr] = fitStateExpr(E1,E2,L)
    fitExpr=L-E2;
end

