function [ nextE_2 ] = StateEqn2( E2,uD1,uD2,L ) % Input: E2(t), uD1(t), uD2(t), L(t)
  global ALPHA_C; global ALPHA_D; global BETA;
  nextE_2=BETA(2)*E2+ALPHA_C(2)*[uD1+uD2-L]-1/ALPHA_D(2)*uD2;
end