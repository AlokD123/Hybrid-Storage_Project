function [ cost ] = CtrlCost_Modified( D1,D2,C2,L ) % Input: D1(t), D2(t), C2(t), L(t)
  global ALPHA_C; global ALPHA_D; global K;
  global MAX_DISCHARGE; global E_MAX; global MAX_CHARGE;
  
  if(D1>MAX_DISCHARGE(1)||D2>MAX_DISCHARGE(2)||D1<0||D2<0||C2>MAX_CHARGE(2)||C2<0) %If control is too high or low
    cost=Inf;                                                                      %Make control inadmissible (infinite cost)
  else                                                                             %Else calculate cost of control
    cost=(1/ALPHA_D(1)-1)*D1/E_MAX(1)+(1/ALPHA_D(2)-1)*D2/E_MAX(2)+(1-ALPHA_C(2))*C2/E_MAX(2)+(1-ALPHA_C(1))*(D1+D2-C2-L)/E_MAX(1)  +K(1)*(D1/E_MAX(1))^2+K(2)*((D1+D2-C2-L)/E_MAX(1))^2;
  end
end