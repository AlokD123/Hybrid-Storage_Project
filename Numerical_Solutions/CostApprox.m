CostStructMtx=[]; %Create matrix containing cost structure
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



%Test selected basis functions by LSQ fitting
%Bases: 1) Constant, 2) E1+E2, 3) L-E2, 4) E2

Phi=[]; %Create design matrix

for i=1:length(CostMtx)
   % Cols contain basis values in each FEASIBLE state (cost~=inf)
   if CostMtx(i,3)<1e9 %For each feasible state...
      phi_vec=[1,E1s(i)+E2s(i),Ls(i)-E2s(i),E2s(i)]; %Create design vector
      Phi=[Phi;phi_vec]; %Add to Phi
   end
end

%Do LSQ on overdetermined system to find r...
r=(Phi'*Phi)^(-1)*Phi'*cost;