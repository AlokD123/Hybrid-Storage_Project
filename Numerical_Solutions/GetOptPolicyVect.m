%PART 1: Create vector of probabilities of states (marginalized over control applied (denominator))
%1) Augment optD vector to include probabilities of infeasible states too (0's)
%Make each E_Ind_Mtx same size to compare ALL states between different
%values of p, once in vector form
%->Augmented vectors for given values of p (i.e. like augmented versions of
%E_Ind_VectALL, and subsets of E_MtxALL_Vect)
E_MtxALL_Vect_subs={};
for p=1:P1*P2
    E_Ind_Mtx_p=E_Ind_Mtx{p};
    E_Ind_Mtx_p(:,size(E_Ind_Mtx_p,2)+1:size(E_Ind_MtxALL,2))=0; %Pad with zeros on side to make same size
    %Convert to vector
    trE_Ind_Mtx_p=E_Ind_Mtx_p';
    E_MtxALL_Vect_subs{p}=trE_Ind_Mtx_p(:);
end

%2) Create augmented vectors of probabilities for ALL states - feasible
%AND INFEASIBLE TOO - for EACH CONTROL p

%For each element in E_MtxALL_Vect_subs{p}, if...
%a) 0, append 0 to aug_optD_subP{p}
%b) non-zero, append some value from optD to aug_optD_subP{p}
%where some value is next value in for i=1:p-1 sumLen=sumLen+len(vecti) end optD(sumLen:sumLen+len(vectp))
aug_optD_subP={};
indOptD=1;
for p=1:P1*P2
    aug_optD_subP_p=[];
    E_MtxALL_Vect_subs_p=E_MtxALL_Vect_subs{p};
   for i=1:length(E_MtxALL_Vect_subs_p)
       if(E_MtxALL_Vect_subs_p(i)==0)
          aug_optD_subP_p=[aug_optD_subP_p;0];
       else
           %Find value in subvector of optD just by continuously
           %indexing through optD in order <--------------------------Assuming optD linearly indexed in order (E2, E1, L, D2, D1)
           aug_optD_subP_p=[aug_optD_subP_p;optD(indOptD)];
           indOptD=indOptD+1;
       end
   end
   aug_optD_subP{p}=aug_optD_subP_p;
end

%3) Marginalise: sum vector components
d_state=zeros(length(aug_optD_subP{1}),1); %Initialize
for(p=1:P1*P2)
   d_state=d_state+aug_optD_subP{p}; %Sum over control values
end

%PART 2: Get stationary probabilities vector
%Create augmented optD vector, for ALL states
aug_optD=[];
for p=1:P1*P2
    aug_optD=[aug_optD;aug_optD_subP{p}];
end
%Create vector with vector d_state duplicated P1*P2 times and appended
%(to allow for dividing each probability by marginalized value)
dup_ones=ones(P1*P2,1);
d_state_dup=kron(dup_ones,d_state);
%Divide to get stationary probabilities vector for ALL states (augmented)
aug_pi=aug_optD./d_state_dup;

%Create augmented vector of all E_MtxALL_Vect_subs vectors
aug_E_MtxALL_Vect=[];
for p=1:P1*P2
    aug_E_MtxALL_Vect=[aug_E_MtxALL_Vect;E_MtxALL_Vect_subs{p}];
end
%Get stationary probabilities vector for ONLY feasible states (non-zero
%in aug_E_MtxALL_Vect)
pi=aug_pi(aug_E_MtxALL_Vect~=0);

%NOTE: in the k-th iteration, policy PI will be the optimal
%policy used later