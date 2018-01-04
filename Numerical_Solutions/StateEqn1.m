function [ nextE_1 ] = StateEqn1( E1,uD1,beta ) % Input: E1(t), uD1(t), and value of beta to use
  global ALPHA_D; global BETA;
  nextE_1=beta*E1-1/ALPHA_D(1)*uD1;
end