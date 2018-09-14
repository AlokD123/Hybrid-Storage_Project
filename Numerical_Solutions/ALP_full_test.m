%Full Approximate LP testing, including interpolation

ApproxLP_sol_IHDP_v12;
GetCtrlPolicy_OptQVals;
PolicyCmpltnessEval_v2;
sum(NumOptCtrls_2(:,4)) %To confirm each state has one single optimal control

LP_IHDP_EvaluatePolicy_v2;
VisualizeIHDP_LP_cost;