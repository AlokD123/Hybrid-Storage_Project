function [ nextE_2 ] = StateEqn2( E2,uD1,uD2,L,beta,alpha_C,alpha_D ) % Input: E2(t), uD1(t), uD2(t), L(t), and value of efficiencies to use
  nextE_2=beta*E2+alpha_C*[uD1+uD2-L]-1/alpha_D*uD2;
end