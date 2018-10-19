function [ nextE_1 ] = StateEqn1( E1,uD1,beta,alpha ) % Input: E1(t), uD1(t), and values of efficiencies to use
  nextE_1=beta*E1-1/alpha*uD1;
end