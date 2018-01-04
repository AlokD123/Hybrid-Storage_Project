function [ nextE_2 ] = StateEqn2( E2,uD1,uD2,L,beta ) % Input: E2(t), uD1(t), uD2(t), L(t), and value of beta to use
  global ALPHA_C; global ALPHA_D;
  nextE_2=beta*E2+ALPHA_C(2)*[uD1+uD2-L]-1/ALPHA_D(2)*uD2;
end