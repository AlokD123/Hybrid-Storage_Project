%Script to TEST BEST APPROXIMATION, for one set of basis functions

global CostMtx;

%% Create matrix containing cost structure

CostStructMtx=[];
E1s=[]; E2s=[]; Ls=[]; %Store values of E1, E2 and L for each state

for i=1:length(CostMtx)
   % Col 1 contains state index E_Ind
   % Cols 2-3 contains associated value of E1 and E2 in each state
   % Cols 4-5 contain load and L-E2
   % Col 6 contains associated cost in each state
   E2=remainder(CostMtx(i,1),N2);
   E1=(CostMtx(i,1)-E2)/N2 + 1;
   L=CostMtx(i,2);
   CostStructMtx(i,:)=[CostMtx(i,1),E1,E2,L,L-E2,CostMtx(i,3)]; %Create matrix
   %Also add values to vectors
   E2s=[E2s;E2]; E1s=[E1s;E1]; Ls=[Ls;L];
end

 %Sort matrix by L-E2
[~,idx] = sort(CostStructMtx(:,5)); % sort just the third column
sortedCostStructMtx = CostStructMtx(idx,:);   % sort the whole matrix using the sort indices



%% Test selected basis functions by 1-norm fitting
%Bases: 1) Constant, 2) E1+E2, 3) L-E2, 4) E2, 5) E1^2, 6) E2^2, 7)
%(L-E2)^2, 8) (L-E1)^2

Phi=[]; %Create design matrix

for i=1:length(CostMtx)
   % Cols contain basis values in each FEASIBLE state (cost~=inf)
   if CostMtx(i,3)<1e9 %For each feasible state...
       %Create parameter fitting vector
      %phi_vec=[1,E1s(i)+E2s(i),Ls(i)-E2s(i),E2s(i)]; 
      %phi_vec=[1,E1s(i)+E2s(i),Ls(i)-E2s(i),E2s(i),E1s(i)^2,E2s(i)^2,(Ls(i)-E2s(i))^2,(Ls(i)-E1s(i))^2, Ls(i)^2,(E2s(i)-E1s(i))^2, Ls(i)^3,E1s(i)^3, E2s(i)^3 ];
      phi_vec=[1,E1s(i)+E2s(i),Ls(i)-E2s(i),E2s(i), E1s(i)^2,E2s(i)^2,(Ls(i)-E2s(i))^2,(Ls(i)-E1s(i))^2, Ls(i)^2,(E2s(i)-E1s(i))^2, Ls(i)^3,E1s(i)^3, E2s(i)^3, E1s(i)^4, E2s(i)^4, Ls(i)^4 ];
      Phi=[Phi;phi_vec]; %Add to Phi
   end
end

%Do 1-norm minimization on overdetermined system to find r...

cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable r_fit(size(Phi,2))
    dual variables d
    minimize( c_state'*(cost-Phi*r_fit) )
    subject to
        d : Q*Phi*r_fit <= b
        Phi*r_fit >= 0
cvx_end

figure

plot(Phi*r_fit,Phi*r_fit, '.');
hold on;
plot(Phi*r_fit,cost, '.');

% optD2=optD;
% optD2(optD2<1e-2)=0;
%  cvx_begin
%     grbControl.LPMETHOD = 1; % Use dual simplex method
%     variable r_fit(size(Phi,2))
%     dual variables d
%     minimize( c_state'*abs(cost-Phi*r_fit) + optD2'*(Q*Phi*r_fit - b) )
% cvx_end

%OLD Do LSQ on overdetermined system to find r...
%r=(Phi'*Phi)^(-1)*Phi'*cost;