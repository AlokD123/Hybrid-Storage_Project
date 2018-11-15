function [formatMtxCosts] = FormatCostVect(cost_vect)
%FormatCostVect Format cost vector into E1xE2 matrices (one for each value of load)
%   Input: vector of costs (actual or approximate)
%   Output: formatted matrices

global E_Ind_MtxALL; global N2; global CostMtx; global MAX_LOAD; global MIN_LOAD;

%1) Create matrix with 3 columns: a) E-state, b) load, c) associated cost
  trE_Ind_MtxALL=E_Ind_MtxALL';
  E_MtxALL_Vect=trE_Ind_MtxALL(:); %Vector form of E_Ind_MtxALL, including infeasible states (0)
  
  CostMtx=[]; %Abovementioned matrix with 3 columns
  indL=0; %Load index, for second column of matrix
  costInd=1; %To index optimal cost sol'n vectors
  indE=0; %To index E-state
  
  for i=1:length(E_MtxALL_Vect)
      %Get E-state value in any case
      if rem(i,MAX_LOAD-MIN_LOAD+1)==1 %New state after cycled through CONSTANT # of loads
          indE=indE+1;
      end
      indL=rem(i,MAX_LOAD-MIN_LOAD+1); %Cycling through constant # of loads
      
      if E_MtxALL_Vect(i)==0       %If infeasible load state...
          %Create row of matrix
          CostMtx(i,:)=[indE,indL,Inf]; %SET INFINITE COST FOR INFEASIBLE STATES
      else                      %OTHERWISE...
          %Create row of matrix
          CostMtx(i,:)=[E_MtxALL_Vect(i),indL,cost_vect(costInd)]; %Get cost from optimal cost solution
          %Take from next row of cost vector if subsequent state feasible
          costInd=costInd+1;
      end
  end
  
  %2) Sort matrix based on loads (i.e. 2nd column)
  sortedCostMtx=sortrows(CostMtx,2);
  
  %3) Get sub-vectors of costs for each distinct load value
  subVectCosts_Load={};
  for i=1:max(sortedCostMtx(:,2)) %For each load value (i)
      subVectCosts_Load{i}=sortedCostMtx(sortedCostMtx(:,2)==i,3); %Get vector of costs with load i
  end
  
  %4) Convert each sub-vector into a set of matrices (E1xE2) called FormatMtxCosts
  formatMtxCosts=[];
  for i=1:length(subVectCosts_Load)
      subVectCosts_Load_i=subVectCosts_Load{i};
      formatMtxCosts(:,:,i)=vec2mat_custom(subVectCosts_Load_i,N2);
  end

end

