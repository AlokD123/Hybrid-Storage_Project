function [ nextE_1 ] = StateEqn1( E1,uD1) % Input: E1(t), uD1(t)
  global ALPHA_D; global BETA;
  nextE_1=BETA(1)*E1-1/ALPHA_D(1)*uD1;
end