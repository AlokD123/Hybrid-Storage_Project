function [ nextE_2 ] = StateEqn2_wRegenBraking_v2( E2,U,L,beta,alpha_C,alpha_D ) % Input: E2(t), U(t), L(t), and value of efficiencies to use
% V2: COMBINED CONTROLS (default mutual exclusion)
  nextE_2=beta*E2-alpha_C*min(0,(L-U))-1/alpha_D*max((L-U),0);
end