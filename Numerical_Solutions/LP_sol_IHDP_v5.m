%LP solution of IHDP (Value Iteration) for Hybrid Storage optimization
%warning('off', 'Octave:possible-matlab-short-circuit-operator');
clearvars -except X V;

global E_MIN; global E_MAX; 
E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
E_MAX=[5;4]; %Maximum energy to be stored (upper bound)

%Input: initial state, horizon
%Initial stored energy (user-defined)
%Must be between MIN_STATE and MAX_STATE
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);

%NOTE: at end, uOpt will have best control policy, and NetCost will contain total cost of DP operation
%ASSUMING FINAL COST = 0


%Model setup
global MAX_CHARGE; global MAX_DISCHARGE;
MAX_CHARGE=[0;100]; %Maximum charging of the supercapacitor
MAX_DISCHARGE=[5;4]; %Maximum discharging of the 1) battery and 2) supercap

global MIN_LOAD;
MIN_LOAD=0; %Minimum load expected
MAX_LOAD=MAX_DISCHARGE(1)+MAX_DISCHARGE(2);

MAX_NUM_ZEROS=3; %Maximum number of zero load counts before end sim

global ALPHA_C; global ALPHA_D; global BETA; global K;
ALPHA_C=[0.99;0.99]; %Efficiency of charging
ALPHA_D=[0.9;0.95]; %Efficiency of discharging
BETA=[0.99;0.99];    %Storage efficiency
K=2;           %Weighting factor for D1^2 cost
PERFECT_EFF=0;

%Discounted infinite horizon problem
global DISCOUNT; %Discount factor
DISCOUNT=[];
DISCOUNT=0.99;


%Definitions
M=MAX_LOAD-MIN_LOAD+1;
N1=(E_MAX(1)-E_MIN(1)+1);
N2=(E_MAX(2)-E_MIN(2)+1);
P1=MAX_DISCHARGE(1)+1;
P2=MAX_DISCHARGE(2)+1;
INF_COST=1000; %Cost of infeasible states (arbitrary sentinel value)

%Initialization
E_Ind_Vect_p=[];      %Vector of current state energies
nextE_Ind_Vect_p=[];  %Vector of next state energies
aug_nextE_Ind_Vect_p=[]; %Augmented vector containing current state energies and next energies for currently infeasible states
numAdmissibleLoads=0; %Count number of admissible load values for a given energy state (for UNIFORM DISTRIBUTION)

PF={};
P_mtx={};
P=[];

indL_Feas=[]; %Vector of feasible demands for ONE GIVEN combination of x and u

unrepNextE_Inds=[]; %List of unrepeated nextE_Ind values

Lmin_p=[]; %Vector of minimum loads required at high discharge (for given p)
Lmin_offs_p=[]; %Vector of minimum load offsets for each E-state, to create CORRECT MAPPING in G matrix

boolDiffNxtState=0; %Flag to indicate different next state different, so don't add current E-state

%% PART A: SET UP MATRICES
%For each possible control...
  for D1=0:MAX_DISCHARGE(1)
    for D2=0:MAX_DISCHARGE(2)
        %Map control to control index
        D1_Ind=D1+1; D2_Ind=D2+1;
        %Get combination #(p)
        p=D2_Ind+P2*(D1_Ind-1);
        
        indCount=0; %Index for feasible state #, for a given value of p
        
        %For each state at an iteration...
        for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
            for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
                %Map state index to state
                E1=E_MIN(1)+(E_Ind1-1);
                E2=E_MIN(2)+(E_Ind2-1);
                
                %Get index of current state energies in vector of state energies
                E_Ind=(E_Ind1-1)*N2+E_Ind2;
                
                if(D1>E1 || D2>E2)  %If discharge too high for state...
                    %IGNORE
                else
                    %Index row in E-state indices mtx (for feasible E-state) same as VALUE of E-state index
                    rowInd_Emtx = E_Ind;
                    %Determine MINIMUM required load for high discharge, to
                    %not overflow E2 (one for each E-state)
                    minL=max(  ceil(1/ALPHA_C(2)*(BETA(2)*E2-E_MAX(2)-D2/ALPHA_D(2))+D1+D2),  0); %Calculate Lmin
                    Lmin_p=[Lmin_p; minL]; %Create vector
                    
                    %For each perturbation at the CURRENT time...
                    for indL=1:(MAX_LOAD-MIN_LOAD+1)
                        %Map index to value of load
                        L=indL+MIN_LOAD-1;

                        %STEP 0
                        %Calculate the state these values of u and w will lead to, even if
                        %impossible...
                        [nextE1,nextE2]=optNextStateLimited(E1,E2,D1,D2,L);
                        if(D1==MAX_DISCHARGE(1)) %<---------------------------------------------------------------------------- SOL#2 for excess discharge: saturate state!!!!!!!!!!!!!!
                           nextE1=0; 
                        end

                        %If next state is amongst those achievable with a given perturbance....
                        if(nextE1<=E_MAX(1) && nextE1>=E_MIN(1))
                            if(nextE2<=E_MAX(2) && nextE2>=E_MIN(2))
                                %IF meeting following conditions: (C_MAX and C_MIN)
                                %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
                                if(~((D1+D2-L)<0||(D1+D2-L)>MAX_CHARGE(2)))
                                  %Count the number of feasible states for a given set of controls (D1,D2)
                                  indCount=indCount+1; %... and use as an index
                                  
                                  %STEP 1
                                  %Add state energy index to vector of FEASIBLE state energies for current value of p (D1,D2 combo)
                                  E_Ind_Vect_p=[E_Ind_Vect_p;E_Ind];
                                  %Add state energy index to matrix of ALL FEASIBLE energies
                                  %DO NOT RESET at end. Will overwrite with same values (and add) each time, which is ok.
                                  E_Ind_MtxALL(rowInd_Emtx,indL)=E_Ind;
                                  E_Ind_Mtx_p(rowInd_Emtx,indL)=E_Ind;
                                    
                                  %Map state to state index, to find cost of next state based on its index
                                  nextE_Ind1=round(nextE1-E_MIN(1)+1);
                                  nextE_Ind2=round(nextE2-E_MIN(2)+1); 

                                  %STEP 2
                                  %Get index of next state energy in vector of state energies
                                  nextE_Ind=(nextE_Ind1-1)*N2+nextE_Ind2;
                                  %Add next state energy index to vector of FEASIBLE next state energies
                                  nextE_Ind_Vect_p=[nextE_Ind_Vect_p;nextE_Ind];
                                  
                                  
                                  
%                                   if(length(nextE_Ind_Vect_p)==1 && E_Ind==nextE_Ind) %If first state being added, and not differing...
%                                       aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %By default, add to augmented vector
%                                       %aug_Lmin_offs_p=[aug_Lmin_offs_p;minL]; %Create vector of minimum load values for each nextE-state, WITH repeats (to add OFFSETS in F matrix)
%                                   else                            %ALL OTHER CASES
%                                       if(above_E_Ind==E_Ind) %If repeating E_Ind...
%                                           %if(boolUpdTemp) 
%                                           unrepNextE_Inds=[unrepNextE_Inds; nextE_Ind]; %Store nextE_Ind value for later (whether currently infeasible OR POSSIBLY NOT)
%                                             %boolUpdTemp=0; %Don't update temp until added this state to aug vector
%                                           %end
%                                           if (boolDiffNxtState==0) %If same next state INITIALLY (i.e. for current E-state with load=0)
%                                             aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %Just add CURRENT E-state to augmented vector
%                                             %aug_Lmin_offs_p=[aug_Lmin_offs_p;minL]; %Add to vector
%                                           else
%                                               %Don't add
%                                           end
%                                       else
%                                           %Else, if adding new E_Ind value (not repeating for load)...
%                                           %For each next state not (/not yet) in aug_nextE_Ind_Vect_p...
%                                           for unrepI=unrepNextE_Inds'
%                                               if ~any(nnz(aug_nextE_Ind_Vect_p==unrepI)) %If NOT in state space of current E-states...
%                                                 for i=1:max(nnz(E_Ind_Vect_p==unrepI),1) %Determine x, number of times to insert (number of loads, at least including load=0)
%                                                     %^--------------ASSUMPTION: E_Ind_Vect_p contains x times, since nextEIndVect increasing in value up to unrepI, and EIndVect(i)>=nextEIndVect(i)
%                                                     %Insert in between (x times)
%                                                    aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p; unrepI];
%                                                    %E2_unrepI=E_MIN(2)+(remainder(unrepI,N2)-1); %Determine E2 energy associate with this E-state
%                                                    %unrep_minL=max(  ceil(1/ALPHA_C(2)*(BETA(2)*E2_unrepI-E_MAX(2)-D2/ALPHA_D(2))+D1+D2),  0); %Calculate Lmin for this unrepeated E-state
%                                                    %aug_Lmin_offs_p=[aug_Lmin_offs_p;unrep_minL]; %Add to vector
%                                                 end
%                                               end
%                                           end
%                                           if (E_Ind==nextE_Ind) %IF should include current E-state (because amongst next E-states)...
%                                             aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %Add regular state index at end, as usual
%                                             %aug_Lmin_offs_p=[aug_Lmin_offs_p;minL]; %Add to vector
%                                             boolDiffNxtState=0;
%                                           else
%                                              boolDiffNxtState=1;
%                                           end
%                                           %boolUpdTemp=1; %Resume updating temp
%                                           unrepNextE_Inds=nextE_Ind; %Restart list of unrepeated nextE_Ind values
%                                       end
%                                   end
%                                   above_E_Ind=E_Ind; %Update "above" E_ind value to current position
                                  
                                  
                                  %STEP
                                  %Add indL to list of FEASIBLE loads for this combination of u and x
                                  indL_Feas=[indL_Feas;indL];
                                  %Create vector of minimum load values for each E-state, WITH repeats (to add OFFSETS in G matrix)
                                  Lmin_offs_p=[Lmin_offs_p;minL];
                                  
                                else
                                  %If no feasible state for this combination of (E1,E2) and L...
                                  nextE_Ind=-1; %Flag next state as impossible
                                end
                            else
                                %If no feasible state for this combination of (E1,E2) and L...
                                nextE_Ind=-1; %Flag next state as impossible
                            end
                        else
                            %If no feasible state for this combination of (E1,E2) and L...
                            nextE_Ind=-1; %Flag next state as impossible
                        end
                        
                        %STEP 4
                        %Create p-th vector g, for constraint
                        if(nextE_Ind~=-1) %If this state leads to a feasible next state...
                            gVec_p(indCount)=CtrlCost(D1,D2,L); %Cost of stage is given by CtrlCost
                        else %Else if infeasible next state...
                            %DO NOTHING
                        end
                    end

                    %Reset feasible loads count, for subsequent energy state
                    numAdmissibleLoads=0;
                    %Reset list of feasible loads (next state)
                    indL_Feas=[];
                end
            end
        end
    
%     %At end of p-th value cycle, if remaining unadded states, add to end
%     for unrepI=unrepNextE_Inds'
%       if ~any(nnz(aug_nextE_Ind_Vect_p==unrepI)) %If NOT in state space of current E-states...
%         for i=1:max(nnz(E_Ind_Vect_p==unrepI),1) %Determine x, number of times to insert (number of loads, at least including load=0)
%             %Insert in between (x times)
%            aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p; unrepI];
%            %E2_unrepI=E_MIN(2)+(remainder(unrepI,N2)-1); %Determine E2 energy associate with this E-state
%            %unrep_minL=max(  ceil(1/ALPHA_C(2)*(BETA(2)*E2_unrepI-E_MAX(2)-D2/ALPHA_D(2))+D1+D2),  0); %Calculate Lmin for this unrepeated E-state
%            %aug_Lmin_offs_p=[aug_Lmin_offs_p;unrep_minL]; %Add to vector
%         end
%       end
%     end
        
    %At end of p-th cycle, restart list of unrepeated nextE_Ind values
    unrepNextE_Inds=[];
    

    %STEP
    %Create augmented vector containing current E-states - EXCLUDING those nextly infeasible - AND ALSO next E-states
    augVectRow=1; %Index row in new augmented vector
    for r=1:length(nextE_Ind_Vect_p) %For each next E-state WITH CURRENT CONTROL COMBO (p)
        %Determine TOTAL number of possible loads for that E-state, given ANY POSSIBLE control used
        numRepNextE=nnz(E_Ind_VectALL==nextE_Ind_Vect_p(r)); %Number of possible loads is number of times repeated in E_Ind_VectALL
        %Add given E-state to augmented vector that many times (for each load)
        aug_nextE_Ind_Vect_p(augVectRow:(augVectRow+numRepNextE-1))=nextE_Ind_Vect_p(r);
    end
    
    %Also, exclude from the augmented vector states that are nextly infeasible
    nextlyInfE=~ismember(aug_nextE_Ind_Vect_p,nextE_Ind_Vect_p);
    aug_nextE_Ind_Vect_p(nextlyInfE)=[];
    Lmin_offs_p(nextlyInfE)=[];
    
    
    
    %Store vector data in cell array
    g{p}=gVec_p';
    E_Ind_Vect{p}=E_Ind_Vect_p;
    nextE_Ind_Vect{p}=nextE_Ind_Vect_p;
    aug_nextE_Ind_Vect{p}=aug_nextE_Ind_Vect_p;
    Lmin{p}=Lmin_p;
    Lmin_offs{p}=Lmin_offs_p;
    %aug_Lmin_offs{p}=aug_Lmin_offs_p;
    
    %Reset matrices/vectors
    nextE_Ind_Vect_p=[];
    E_Ind_Vect_p=[];
    aug_nextE_Ind_Vect_p=[];
    gVec_p=[];
    Lmin_p=[];
    Lmin_offs_p=[];
    %aug_Lmin_offs_p=[];
    
    numAdmissibleLoads=0;
    
    end
  end
  
  
  %STEP : Construct vector of ALL FEASIBLE energies
  E_Ind_VectALL=[];
  for row=1:size(E_Ind_MtxALL,1)
      nnzRow=nnz(E_Ind_MtxALL(row,:));
      E_Ind_Mtx_nzRow=E_Ind_MtxALL(row,1:nnzRow);
      E_Ind_VectALL=[E_Ind_VectALL; E_Ind_Mtx_nzRow'];
  end
  
  %STEP 
  %Create full probability matrix
  %SET DISTRIBUTION: UNIFORM
  %(Note: can't create until E_Ind_MtxALL complete, so outside main loop)
  for r=1:size(E_Ind_MtxALL,1)
      P_fullmtx(r,:)=E_Ind_MtxALL(r,:)/sum(E_Ind_MtxALL(r,:)); %<----------------For UNIFORM probability, just NORMALIZE rows of feasible states!!
  end
  
  %STEP 6
  %Create P matrix: select rows corresponding to components in nextE_Ind_Vect
  %(Doing after P_fullmtx completed)
  for p=1:P1*P2
    E_Ind_Vect_p=E_Ind_Vect{p};
    nextE_Ind_Vect_p=nextE_Ind_Vect{p};
    aug_nextE_Ind_Vect_p=aug_nextE_Ind_Vect{p};
    
    for r=1:length(E_Ind_Vect_p)
        Ind_nextE=nextE_Ind_Vect_p(r);    %Get index of state stored in r-th row of nextE_Ind_Vect (i.e. the next energy state)
        
        %Get column number of next row of probabilities as RELATED to the NEXT ENERGY STATE INDEX (mapping to deterministic component!!!)
        c=find(aug_nextE_Ind_Vect_p==Ind_nextE,1); %Get from position of FIRST Ind_nextE in AUG_nextE_Ind_Vect!!!!! (b/c same width as AUGMENTED VECTOR)
        
        %Count number of non-zero probabilities in associated E-state row of P_fullmtx (i.e. Ind_nextE)
        nnzProb_nextE=nnz(P_fullmtx(Ind_nextE,:));      %Should be equal to number of repeats in nextE_Ind_Vect
        %Get said non-zero probabilities
        prob_nextE=nonzeros(P_fullmtx(Ind_nextE,:));
        
        %Fill in row r with said probabilities
        P(r,c:(c+nnzProb_nextE-1))=prob_nextE';
    end
        
    %Store in p-th PF matrix, as well as in own P_mtx
    PF{end+1}=P; 
    %RESET
    P=[];
  end
  
  %STEP : Construct each F matrix
  for p=1:P1*P2
      aug_nextE_Ind_Vect_p=aug_nextE_Ind_Vect{p};
      %aug_Lmin_offs_p=aug_Lmin_offs{p};
      %Index COLUMN of F matrix by ROW number of E_Ind_VectALL
      row=1; %Reset row being checked in E_Ind_VectALL to start when start on next E_Ind vector
      
      %Go through next E-state index vector for current value of p...
      for r=1:length(aug_nextE_Ind_Vect_p)
          %If next state is currently infeasible...
          if aug_nextE_Ind_Vect_p(r)<E_Ind_VectALL(row) %(i.e. NOT continuously increasing in augmented vector)
             row=1; %Restart from beginning of E_Ind_VectALL to find the state <----- ASSUMING ONLY 1 distinct new currently infeasible state!
          end
          
          %Find distinct new nextE-state
          if(r==1)
              boolNewNextEState=1;
          else
              if(aug_nextE_Ind_Vect_p(r)~=aug_nextE_Ind_Vect_p(r-1))
                 boolNewNextEState=1;
              else
                  boolNewNextEState=0;
              end
          end
          
          while(E_Ind_VectALL(row)~=aug_nextE_Ind_Vect_p(r)) %While not reached mapping column in F (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          
          F_p(r,row)=1; %Once reached, map
          row=row+1;  %Start at next column in F next time <-------- Assuming continuously increasing in augmented vector (fixed above)
      end
      
      if isempty(F_p)   %If empty, ignore
         F_p=0;
      else      %IN MOST CASES...
        %Add extra zeros at end to ensure dimensions of F_p and E_Ind_VectALL match
        F_p(:,(size(F_p,2)+1):length(E_Ind_VectALL))=0;
      end
      
      F{p}=F_p;
      F_p=[]; %Reset
      
      if isempty(PF{p}) %If no next state..
          PF{p}=0;  %Ignore constraint
      else      %IN MOST CASES...
         PF{p} = PF{p}*F{p}; %Finish PF matrices 
      end
  end
  
  %STEP : Construct each G matrix
  for p=1:P1*P2
      E_Ind_Vect_p=E_Ind_Vect{p};
      Lmin_offs_p=Lmin_offs{p};
      %Index COLUMN of G matrix by ROW number of E_Ind_VectALL
      row=1; %Reset row being checked in E_Ind_VectALL to start when start on next E_Ind vector
      
      %Go through E-state index vector for current value of p...
      for r=1:length(E_Ind_Vect_p)
          %Find distinct new E-state
          if(r==1)
              boolNewEState=1;
          else
              if(E_Ind_Vect_p(r)~=E_Ind_Vect_p(r-1))
                 boolNewEState=1;
              else
                  boolNewEState=0;
              end
          end
          while(E_Ind_VectALL(row)~=E_Ind_Vect_p(r)) %While not reached mapping column in G (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          if(boolNewEState==1)  %Only if distinct new state...
              row=row+Lmin_offs_p(r); %Add minimum load offset to first state #
          else
              %Otherwise, do nothing because already starting from offset
          end
          G_p(r,row)=1; %Once reached, map
          row=row+1;  %Start at next column in G next time (continuously increasing)
      end
      
      if isempty(G_p)   %If empty, ignore constraint
         G_p=0;
      else      %IN MOST CASES...
        %Add extra zeros at end to ensure dimensions of G_p and E_Ind_VectALL match
        G_p(:,(size(G_p,2)+1):length(E_Ind_VectALL))=0;
      end
      
      G{p}=G_p;
      G_p=[]; %Reset
  end
  
  %Set g=zeros for empty g vectors
  for p=1:P1*P2
    if isempty(g{p})
       g{p}=zeros(length(E_Ind_VectALL),1);
    end
  end
  
  %% CORRECTED MATRICES (REMAINING infeasible states removed)
  %1) Coefficients
  %Create full 'A' matrices for coefficients (A=G-alpha*PF)
  for p=1:P1*P2
    A{p}=G{p}-DISCOUNT*PF{p};
  end
  %Adjoin A matrices to form Q
  Q=[];
  for p=1:P1*P2
    Q=[Q;A{p}];
  end
  %If empty columns in Q...
  if (~all(any(Q,1)))
      disp('ERROR!!!!!'); %ERROR
     Q(:,~any(Q,1))=[]; %Remove, for now 
  end
  
  %2) Constants
  %Create full 'b' vector for constants
  b=[];
  for p=1:P1*P2
    b=[b;g{p}];
  end
  
  %% PART B: OPTIMIZATION
  %Created LP matrices and vectors.
  %Run optimization problem, and find primal as well as dual.
  cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable cost(size(Q,2))
    dual variables d
    maximize( sum(cost) )
    subject to
        d : Q*cost <= b
  cvx_end
  %Get vector of optimal dual
  optD = d;
  
  %Format cost vector into E1xE2 matrices (M matrices, for each value of load)
    for i1=0:M*N2:M*N2*(N1-1)
        for j=0:M-1
          for ind=(1+i1+j):M:(M*(N2-1)+i1+(j+1))
              if(mod(ind,M*N2)==0)
                ind2=(M*N2-1-j)/M+1;
              else
                ind2=(mod(ind,M*N2)-1-j)/M+1;
              end
              ConvCosts(i1/(M*N2)+1,ind2,j+1)=cost(ind);
          end
        end
    end
    
    %% PART C: STATIONARY POLICY
    %Create vector of probabilities marginalized over control applied (denominator)
    d_state=zeros(N1*N2*M,1);
    for(iij=1:N1*N2*M)
       for(p=1:P1*P2)
           d_state(iij)=d_state(iij)+optD(iij+N1*N2*M*(p-1)); %Sum over control values
       end 
    end
    %Create vector with vector d_state duplicated P1*P2 times and appended
    %(to allow for dividing each probability by marginalized value)
    dup_ones=ones(P1*P2,1);
    d_state_dup=kron(dup_ones,d_state);
    %Divide to get stationary probabilities vector
    pi=optD./d_state_dup;
    
    
  %Quantify difference between costs and expected matrix of IHDP costs
  Diff=X-ConvCosts; %Matrix of differences. X comes from Value Iteration
  Diff(Diff==Inf)=0; %Ignore infinite differences (correct)
  %Use 2-norm to quantify difference
  norm_array=arrayfun(@(idx) norm(Diff(:,:,idx),1), 1:size(Diff,3));
  TotalDiff=norm(norm_array,1);
  disp('Net difference is'); TotalDiff