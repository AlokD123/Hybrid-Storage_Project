function [ cost ] = CtrlCost_Modified_v3( U1,L,E1 ) % Input: U1(t), L(t), E1(t)
%V3: UPDATED stage cost
  global ALPHA_C; global ALPHA_D; global K;
  global MAX_DISCHARGE; global E_MAX; global MAX_CHARGE;
  
  global nu; nu=1; %Cost for high state of charge
  global sigma; sigma=1; %To tune cost of low state of charge
  
  if(U1>MAX_DISCHARGE(1)||(L-U1)<-MAX_CHARGE(2)||U1<-MAX_CHARGE(1)||(L-U1)>MAX_DISCHARGE(2)) %If control is too high or low
    cost=Inf;                                                                                %Make control inadmissible (infinite cost)
  else                                                                                       %Else calculate cost of control
    cost=(1/ALPHA_D(1)-1)*max(U1,0)/E_MAX(1)+(1/ALPHA_D(2)-1)*max(L-U1,0)/E_MAX(2)-(1-ALPHA_C(2))*min(0,L-U1)/E_MAX(2)-(1-ALPHA_C(1))*min(0,U1)/E_MAX(1)  +(E_MAX(1)/(E1+sigma))*(K(1)*(max(U1,0)/E_MAX(1))^2+K(2)*(min(0,U1)/E_MAX(1))^2) +nu*E1/E_MAX(1);
  end
end