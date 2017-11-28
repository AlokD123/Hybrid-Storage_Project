function [ nextE ] = StateEqn( E,uD1,uD2,L ) % Input: E(:,t), uD1(t), uD2(t), L(t)
  global ALPHA_C; global ALPHA_D; global BETA;
  nextE_1=BETA(1)*E(1)-1/ALPHA_D(1)*uD1;
  nextE_2=BETA(2)*E(2)+ALPHA_C(2)*[uD1+uD2-L]-1/ALPHA_D(2)*uD2;
  nextE = [nextE_1;nextE_2];
end