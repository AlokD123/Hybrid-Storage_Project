function [ nextE_2 ] = StateEqn2_wRegenBraking( E2,uD2,uC2,beta ) % Input: E2(t), uD2(t), uD2(t), and value of beta to use
  global ALPHA_C; global ALPHA_D;
  nextE_2=beta*E2+ALPHA_C(2)*uC2-1/ALPHA_D(2)*uD2;
end