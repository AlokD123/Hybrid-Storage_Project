function [ cost ] = CtrlCost( D1,D2,L ) % Input: D1(t), D2(t), L(t)
  global ALPHA_C; global ALPHA_D; global BETA; global K;
  cost=(1-ALPHA_D(1))*D1+K*(D1)^2+(1-ALPHA_D(2))*D2+(1-ALPHA_C(2))*(D1+D2-L);
end