%Create matrix of optimal control policy, fullPolicyMtx 
%fullPolicyMtx(E1,E2,L,:) is control vector of state (E1,E2,L), containing 1 in optimal p-th element
%policyMtx contains policy in elements, but also indicates infeasible states (NaN)

clear policyMtx;
global fullPolicyMtx;

%N1=6; N2=5; M=10;
%P1=N1; P2=N2;

%Repeat for each decision (D1,D2 combination)
for(p=1:p_max)
piSubVec=aug_pi((p-1)*N1*N2*(M-1)+1:p*N1*N2*(M-1));
%Format binary policy sub-vector (for cases when applying the single decision) into E1xE2 matrices ((M-1) matrices, for each value of load)
    for(i1=0:(M-1)*N2:(M-1)*N2*(N1-1))
        for(j=0:(M-1)-1)
          for(ind=(1+i1+j):(M-1):((M-1)*(N2-1)+i1+(j+1)))
              if(mod(ind,(M-1)*N2)==0)
                ind2=((M-1)*N2-1-j)/(M-1)+1;
              else
                ind2=(mod(ind,(M-1)*N2)-1-j)/(M-1)+1;
              end
              policyMtx(i1/((M-1)*N2)+1,ind2,j+1,p)=piSubVec(ind);
          end
        end
    end
end

%Create fullPolicyMtx by replacing NaNs with 0s
%Note: NaN states not feasible, so can ignore (0)
fullPolicyMtx=policyMtx;
fullPolicyMtx(isnan(fullPolicyMtx))=0;