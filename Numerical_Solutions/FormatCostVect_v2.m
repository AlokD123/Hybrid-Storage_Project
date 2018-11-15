function [formatMtxCosts] = FormatCostVect_v2(cost_vect)
%FormatCostVect_v2 Format cost vector into E1xE2 matrices (one for each value of load)
%   Input: vector of costs (actual or approximate)
%   Output: formatted matrices

global E_Ind_MtxALL; global N2; global CostMtx;

%1) Create matrix with 3 columns: a) E-state, b) load, c) associated cost
  trE_Ind_MtxALL=E_Ind_MtxALL';
  E_MtxALL_Vect=trE_Ind_MtxALL(:); %Vector form of E_Ind_MtxALL, including infeasible states (0)
  
  CostMtx=[]; %Abovementioned matrix with 3 columns
  indL_minL=0; %Load index, for second column of matrix. THIS IS RELATIVE TO min(MINLoad_E_state), NOT RELATIVE TO MINLoad_E_state for the current E-state!!!!!!
  costInd=1; %To index optimal cost sol'n vectors
  lastNz=1; %Last feasible E-state (row) index
  
  for i=1:length(E_MtxALL_Vect)
      if rem(i,size(E_Ind_MtxALL,2))==1 %If starting on new E-state (multiple of total # of all loads + 1), reset
          indL_minL=0; %Reset to 0
          lastNz=lastNz+1; %Next E-state, even if not where last non-zero
      end
      if E_MtxALL_Vect(i)==0       %If infeasible load state...
          indL_minL=indL_minL+1;      %Next load
          %Create row of matrix
          CostMtx(i,:)=[lastNz,indL_minL,Inf]; %SET INFINITE COST FOR INFEASIBLE STATES
      else                      %OTHERWISE...
          %Determine whether going onto next row of E_Ind_MtxALL
          if i==1 || E_MtxALL_Vect(i)==E_MtxALL_Vect(i-1) %For same E-state (row)
              indL_minL=indL_minL+1;      %Next load
          else
              indL_minL=indL_minL+1; %DO NOT Restart load indexing, even for new E-state
          end
          
          lastNz=E_MtxALL_Vect(i); %Distinct E-state value (row of E_Ind_MtxALL)
          %Create row of matrix
          CostMtx(i,:)=[E_MtxALL_Vect(i),indL_minL,cost_vect(costInd)]; %Get cost from optimal cost solution
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

