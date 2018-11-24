function [D1,D2,C2] = GetPOpt_wo_Interp_wRegenB_v2(indE1,indE2,L)
%GetPOpt_wo_Interp_wRegenB_v2 Returns the optimal controls in the given state.
%                             Accounts for R.B.
%                             Accounting for D2,C2 MUTUAL EXCLUSION constraint
%   Input: state, assumed ON GRID
%   Output: optimal control

global fullPolicyMtx; global P2; global P3; global MIN_LOAD; global epsilon3; global epsilon4;
global qValsMtx;

global RES_D2; global RES_C2;

%RES_D2=1; RES_C2=1;

pOpt=find(fullPolicyMtx(indE1,indE2,L-MIN_LOAD+1,:)==1); %(Find equal to 1)

%Calculate the net energy change in the supercapacitor
netE2_ctrl=rem(pOpt,P2+P3-1); %:=C2_Ind+D2_Ind-1

%Get associated optimal controls when accounting for piecewise C2,D2 (mutual exclusion)
if netE2_ctrl<=P3
    C2=(netE2_ctrl-1)/RES_C2;
    D2=0;
else
    C2=0;
    D2=rem(netE2_ctrl,P3)/RES_D2;
end

D1=(pOpt-netE2_ctrl)/(P2+P3-1)+1;

end

