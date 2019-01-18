function U1 = GetPOpt_wo_Interp_wRegenB_v3(indE1,indE2,L)
%GetPOpt_wo_Interp_wRegenB Returns the optimal control in the given state.
%                          Accounts for R.B.
%   Input: state, assumed ON GRID
%   Output: optimal control

global fullPolicyMtx; global MIN_LOAD; global RES_U1; global MAX_CHARGE;
global RES_STATE;

pOpt=find(fullPolicyMtx(indE1,indE2,RES_STATE*(L-MIN_LOAD)+1,:)==1); %(Find equal to 1)

%Get associated optimal control
U1_Ind=pOpt;
U1=(U1_Ind-1)/RES_U1-MAX_CHARGE(1);
end

