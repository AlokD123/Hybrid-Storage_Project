cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable cost(length(PFId))
    dual variables d{P1*P2}
    minimize( -1*sum(cost) )
    subject to
        for p=1:P1*P2
            d{p} : (eye(length(PFId))-DISCOUNT*PFId(:,:,p))*cost <= g(:,p)
        end
  cvx_end