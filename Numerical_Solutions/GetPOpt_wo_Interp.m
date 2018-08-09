function [D1,D2] = GetPOpt_wo_Interp(indE1,indE2,L)
%GetPOpt_wo_Interp Returns the optimal controls in the given state
%   Input: state, assumed ON GRID
%   Output: optimal control

global fullPolicyMtx; global P2; global MIN_LOAD;

pOpt=find(abs(fullPolicyMtx(indE1,indE2,L-MIN_LOAD+1,:)-1)<0.01); %(Find equal to 1)

%Get associated optimal controls
D2_Ind=remainder(pOpt,P2);
D1_Ind=(pOpt-D2_Ind)/P2+1;
D1=D1_Ind-1; D2=D2_Ind-1;

end

