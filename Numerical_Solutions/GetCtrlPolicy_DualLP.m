%Create matrix of optimal control policy, fullPolicyMtx 
%fullPolicyMtx(E1,E2,L,:) is control vector of state (E1,E2,L), containing 1 in optimal p-th element
%policyMtx contains policy in elements, but also indicates infeasible states (NaN)
%For NO regenerative braking case (uncombined controls)

clear policyMtx;
global fullPolicyMtx;

%N1=6; N2=5; M=10;
%P1=N1; P2=N2;

Sigma=size(E_Ind_MtxALL,2); %i.e. -E_MAX -> M-2

%Repeat for each decision (D1,D2 combination)
for p=1:(p_max) 
    if (~isempty(Lmin{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...
        piSubVec=aug_pi((p-1)*N1*N2*Sigma+1:p*N1*N2*Sigma);
        %Format binary policy sub-vector (for cases when applying the single decision) into E1xE2 matrices (Sigma matrices, for each value of load)
            for i=0:Sigma*N2:Sigma*N2*(N1-1) 
                for j=0:Sigma-1 
                  for ind=(1+i+j):Sigma:(Sigma*(N2-1)+i+(j+1)) 
                      if(mod(ind,Sigma*N2)==0)
                        ind2=(Sigma*N2-1-j)/Sigma+1;
                      else
                        ind2=(mod(ind,Sigma*N2)-1-j)/Sigma+1;
                      end
                      policyMtx(i/(Sigma*N2)+1,ind2,j+1,p)=piSubVec(ind);
                  end
                end
            end
    end
end

%Create fullPolicyMtx by replacing NaNs with 0s
%Note: NaN states not feasible, so can ignore (0)
fullPolicyMtx=policyMtx;
fullPolicyMtx(isnan(fullPolicyMtx))=0;