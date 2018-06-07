%Script to TEST BEST APPROXIMATION, for one set of basis functions
%V2: LARGE NUMBER OF BASIS FUNCTIONS (high dimensional model, dimension n)

global CostMtx;
n=8; %Order of polynomial approximation

%% Create matrix containing cost structure

CostStructMtx=[];
E1s=[]; E2s=[]; Ls=[];              %Store values of E1, E2 and L for each state
feasE1s=[]; feasE2s=[]; feasLs=[];  %Store FEASIBLE values of E1, E2 and L for each state

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
   %1) All values
   E2s=[E2s;E2]; E1s=[E1s;E1]; Ls=[Ls;L];
   %2) Feasible values
   if CostMtx(i,3)<1e9 %For each feasible state...
      feasE2s=[feasE2s;E2]; feasE1s=[feasE1s;E1]; feasLs=[feasLs;L]; %Add
   end
end

 %Sort matrix by L-E2
[~,idx] = sort(CostStructMtx(:,5)); % sort just the third column
sortedCostStructMtx = CostStructMtx(idx,:);   % sort the whole matrix using the sort indices



%% Test selected basis functions by 1-norm fitting
%Bases: all terms in order-n polynomial

Phi=[]; %Design matrix to create

%Adjoin feasible state vectors to form a ?x3 array
feasStatesArr=[feasE1s,feasE2s,feasLs];

%Create design matrix with fitting functions up to order n
Phi=DesignMtx(feasStatesArr,cost,n);

%Do 1-norm minimization on overdetermined system to find r...

cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable r_fit(size(Phi,2))
    dual variables d
    minimize( c_state'*abs(cost-Phi*r_fit) )
    subject to
        d : Q*Phi*r_fit <= b
cvx_end

figure
hold on;
plot(Phi*r_fit,Phi*r_fit, '*');
plot(Phi*r_fit,cost, '*');
xlabel('Approximate Cost');
legend('Approximate Cost','Actual Cost');
title(strcat('Testing Cost Approximation Using Order-',num2str(n),' Fit (',num2str(length(r_fit)),' Bases)'));
%plot(Phi*r_fit, '.');
%plot(cost, '.');

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