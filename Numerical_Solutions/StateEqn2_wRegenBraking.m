function [ nextE_2 ] = StateEqn2_wRegenBraking( E2,uD2,uC2,beta,alpha_C,alpha_D ) % Input: E2(t), uD2(t), uD2(t), and value of efficiencies to use
  nextE_2=beta*E2+alpha_C*uC2-1/alpha_D*uD2;
end