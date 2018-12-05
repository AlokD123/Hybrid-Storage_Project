function [ nextE_1 ] = StateEqn1_wRegenBraking_v2( E1,U,beta ) % Input: E1(t), U(t), and value of beta to use
% V2: COMBINED CONTROLS (default mutual exclusion)
  global ALPHA_D; global ALPHA_C;
  nextE_1=beta*E1-ALPHA_C(1)*min(0,U)-1/ALPHA_D(1)*max(U,0);
end