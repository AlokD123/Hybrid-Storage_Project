function [ cost ] = CtrlCost( D1,D2,L ) % Input: D1(t), D2(t), L(t)
  %Calculate staget cost
  %For NO regenerative braking case (uncombined controls)
  global ALPHA_C; global ALPHA_D; global K;
  global MAX_DISCHARGE;
  
  if(D1>MAX_DISCHARGE(1)||D2>MAX_DISCHARGE(2)||D1<0||D2<0) %||(D1+D2)>L %If control is too high or low
    cost=Inf;                                              %Make control inadmissible (infinite cost)
  else
    cost=K*(D1)^2+((1/ALPHA_D(1)-1)*D1+(1/ALPHA_D(2)-1)*D2+(1-ALPHA_C(2))*(D1+D2-L)); %Else calculate cost of control
  end
end