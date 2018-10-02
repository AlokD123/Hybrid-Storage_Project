function [ nextE_1 ] = StateEqn1_wRegenBraking( E1,uD1,uD2,uC2,L,beta ) % Input: E1(t), uD1(t), uD2(t), uC2(t), L(t), and value of beta to use
  global ALPHA_D; global ALPHA_C;
  nextE_1=beta*E1+ALPHA_C(1)*(uD1+uD2-uC2-L)-1/ALPHA_D(1)*uD1;
end