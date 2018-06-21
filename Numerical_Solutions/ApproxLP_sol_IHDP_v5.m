%LP solution of IHDP (Value Iteration) for Hybrid Storage optimization
%USING COST APPROXIMATION
% v5: AUTOMATIC BASIS FUNCTION GENERATION (iterative), v3

%warning('off', 'Octave:possible-matlab-short-circuit-operator');
clearvars -except X V cost optD_exact optD_approx_old;

global E_MIN; global E_MAX;
E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
E_MAX=[5;4]; %Maximum energy to be stored (upper bound)

%Solver tolerance
tolerance=1e-6;

%% Input: initial state, horizon, number of bases in approximation
%Initial stored energy (user-defined)
%Must be between MIN_STATE and MAX_STATE
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);

R=5; %MAXIMUM number of bases (MUST BE LESS THAN NUMBER OF FEASIBLE STATES)

%% Model setup
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


%% Definitions
global N2;

M=MAX_LOAD-MIN_LOAD+1;
N1=(E_MAX(1)-E_MIN(1)+1);
N2=(E_MAX(2)-E_MIN(2)+1);
P1=MAX_DISCHARGE(1)+1;
P2=MAX_DISCHARGE(2)+1;
INF_COST=1000; %Cost of infeasible states (arbitrary sentinel value)

%% Initialization
E_Ind_Vect_p=[];      %Vector of current state energies
nextE_Ind_Vect_p=[];  %Vector of next state energies
aug_nextE_Ind_Vect_p=[]; %Augmented vector containing current state energies and next energies for currently infeasible states
numAdmissibleLoads=0; %Count number of admissible load values for a given energy state (for UNIFORM DISTRIBUTION)

P_mtx={};   %Array of P matrices
P=[];       %Current P matrix
PF={};      %Array of P*F matrices

global E_Ind_MtxALL; %Matrix of all states (0 if infeasible)
global CostMtx; %Matrix with 3 columns: a) E-state, b) load, c) associated cost. Used for approximation
CostMtx=[];

indL_Feas=[]; %Vector of feasible demands for ONE GIVEN combination of x and u
feasStates=[]; %List of all feasible states (E1,E2,L), no repeats

Lmin_p=[]; %Vector of minimum loads required at high discharge (for given p)
Lmin_offs_p=[]; %Vector of minimum load offsets for each E-state, to create CORRECT MAPPING in G matrix
E_Ind_Mtx_p=[]; %Matrix of E_Ind_MtxALL values, but for EACH value of p


Phi=[]; %Design matrix, for cost approximation

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
                                %IF meeting following conditions: (C_MIN and C_MAX)
                                %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
                                if(~((D1+D2-L)<0||(D1+D2-L)>MAX_CHARGE(2)))
                                  %Count the number of feasible states for a given set of controls (D1,D2)
                                  indCount=indCount+1; %... and use as an index
                                  
                                  %STEP 1: create vector and matrix of FEASIBLE state energies for each load
                                  %Add state energy index to vector for current value of p (D1,D2 combo)
                                  E_Ind_Vect_p=[E_Ind_Vect_p;E_Ind];
                                  %Add state energy index to matrix of ALL FEASIBLE energies
                                  %DO NOT RESET at end. Will overwrite with same values (and add) each time, which is ok.
                                  E_Ind_MtxALL(rowInd_Emtx,indL)=E_Ind;
                                  E_Ind_Mtx_p(rowInd_Emtx,indL)=E_Ind;
                                    
                                  %Map state to state index, to find cost of next state based on its index
                                  nextE_Ind1=round(nextE1-E_MIN(1)+1);
                                  nextE_Ind2=round(nextE2-E_MIN(2)+1); 

                                  %STEP 2: create vector of next state energies for each load
                                  %Get index of next state energy in vector of state energies
                                  nextE_Ind=(nextE_Ind1-1)*N2+nextE_Ind2;
                                  %Add next state energy index to vector of FEASIBLE next state energies
                                  nextE_Ind_Vect_p=[nextE_Ind_Vect_p;nextE_Ind];
                                   
                                  %STEP 3: determine feasible loads
                                  %Add indL to list of FEASIBLE loads for this combination of u and x
                                  indL_Feas=[indL_Feas;indL];
                                  %Create vector of minimum load values for each E-state, WITH repeats (to add OFFSETS in G matrix)
                                  Lmin_offs_p=[Lmin_offs_p;minL];
                                  
                                  %STEP 4: Create list of all FEASIBLE states
                                  feasStates(E_Ind1,E_Ind2,indL)=1;
                                  
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
                        
                        %STEP 5
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
        
    %Store vector data in cell array
    g{p}=gVec_p';
    E_Ind_Vect{p}=E_Ind_Vect_p;
    nextE_Ind_Vect{p}=nextE_Ind_Vect_p;
    Lmin{p}=Lmin_p;
    Lmin_offs{p}=Lmin_offs_p;
    E_Ind_Mtx{p}=E_Ind_Mtx_p;
    
    %Reset matrices/vectors
    nextE_Ind_Vect_p=[];
    E_Ind_Vect_p=[];
    gVec_p=[];
    Lmin_p=[];
    Lmin_offs_p=[];
    E_Ind_Mtx_p=[];
    
    numAdmissibleLoads=0;
    
    end
  end
  
  
  %STEP 6: Construct vector of ALL FEASIBLE energies, for all control
  E_Ind_VectALL=[];
  for row=1:size(E_Ind_MtxALL,1)
      nnzRow=nnz(E_Ind_MtxALL(row,:));
      E_Ind_Mtx_nzRow=E_Ind_MtxALL(row,1:nnzRow);
      E_Ind_VectALL=[E_Ind_VectALL; E_Ind_Mtx_nzRow'];
  end
  
  %STEP 7: Create full probability matrix
  %SET DISTRIBUTION: UNIFORM
  %(Note: can't create until E_Ind_MtxALL complete, so outside main loop)
  for r=1:size(E_Ind_MtxALL,1)
      P_fullmtx(r,:)=E_Ind_MtxALL(r,:)/sum(E_Ind_MtxALL(r,:)); %<----------------For UNIFORM probability, just NORMALIZE rows of feasible states!!
  end
  
  
  
  for p=1:P1*P2
    E_Ind_Vect_p=E_Ind_Vect{p};
    nextE_Ind_Vect_p=nextE_Ind_Vect{p};
    Lmin_offs_p=Lmin_offs{p};
    
    %STEP 8: Create augmented vector containing current E-states - EXCLUDING those nextly infeasible - AND ALSO next E-states
    %(Note: doing after E_Ind_VectALL complete)
    augVectRow=1; %Index row in new augmented vector
    r=1; %Start from beginning
    while r<(length(nextE_Ind_Vect_p)+1) %For each next E-state WITH CURRENT CONTROL COMBO (p)
        if (r~=1) %...IN MOST CASES
            %If next E-state already counted once, do not double-count...
            while r<(length(nextE_Ind_Vect_p)+1) && nnz(nextE_Ind_Vect_p(1:r-1)==nextE_Ind_Vect_p(r))
                r=r+1; %Skip to next unrepeated E-state
            end
        end
        if(r~=(length(nextE_Ind_Vect_p)+1))
            %Determine TOTAL number of possible loads for that E-state, given ANY POSSIBLE control used
            numRepNextE=nnz(E_Ind_VectALL==nextE_Ind_Vect_p(r)); %Number of possible loads is number of times repeated in E_Ind_VectALL
            %Add given E-state to augmented vector that many times (for each load)
            aug_nextE_Ind_Vect_p(augVectRow:(augVectRow+numRepNextE-1),1)=nextE_Ind_Vect_p(r);
            augVectRow=augVectRow+numRepNextE; %Start adding at end next time 
        end
        r=r+1; %Manually increment index in while loop
    end
    
    %Also, exclude from the augmented vector states that are nextly infeasible
    nextlyInfE=~ismember(aug_nextE_Ind_Vect_p,nextE_Ind_Vect_p);
    aug_nextE_Ind_Vect_p(nextlyInfE)=[];
    Lmin_offs_p(nextlyInfE)=[];
    
    %Store in cell array
    aug_nextE_Ind_Vect{p}=aug_nextE_Ind_Vect_p;
    
    %STEP 9: Create each P matrix
    %For P matrix, select rows corresponding to components in nextE_Ind_Vect
    %(Note: doing after P_fullmtx completed)
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
    PF{p}=P;
    P_mtx{p}=P;
    %Reset matrices/vectors
    P=[];
    aug_nextE_Ind_Vect_p=[];
  end
  
  
  
  
  %STEP 10: Construct each F matrix
  for p=1:P1*P2
      aug_nextE_Ind_Vect_p=aug_nextE_Ind_Vect{p};
      %Index COLUMN of F matrix by ROW number of E_Ind_VectALL
      row=1; %Reset row being checked in E_Ind_VectALL to start when start on next E_Ind vector
      
      %Go through next E-state index vector for current value of p...
      for r=1:length(aug_nextE_Ind_Vect_p)
          %If next state is currently infeasible...
          if aug_nextE_Ind_Vect_p(r)<E_Ind_VectALL(row) %(i.e. NOT continuously increasing in augmented vector)
             row=1; %Restart from beginning of E_Ind_VectALL to find the state <----- ASSUMING ONLY 1 distinct new currently infeasible state!
          end
          
          while(E_Ind_VectALL(row)~=aug_nextE_Ind_Vect_p(r)) %While not reached mapping column in F (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          
          F_p(r,row)=1; %Once reached, map
          row=min(row+1,length(E_Ind_VectALL)); %Start at next column in F next time, saturating at maximum
          %^-------- Assuming continuously increasing in augmented vector (fixed above)
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
  
  %STEP 11: Construct each G matrix
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
          row=min(row+1,length(E_Ind_VectALL));  %Start at next column in G next time (continuously increasing)
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
  
  %If get empty g vector...
  for p=1:P1*P2
    if isempty(g{p})
        disp('ERROR!!'); %Error!!!
       g{p}=zeros(length(E_Ind_VectALL),1); %Set equal to zeros
    end
  end
  
  %% CORRECTED MATRICES (REMAINING infeasible states removed)
  %1) Coefficients
  Q=[];
  for p=1:P1*P2
    %Create full 'A' matrices for coefficients (A=G-alpha*PF)
    A{p}=G{p}-DISCOUNT*PF{p};
    %Adjoin A matrices to form Q
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
  
  %% PART B: COST APPROXIMATION AND BASIS FUNCTION GENERATION
  %Create initial design matrix (1 row per feasible state)
  Phi=ones(size(Q,2),1);

% Find state-relevance vector for minimization, c
% TAKE c TO BE STEADY STATE ENTERING PROBABILITIES FOR EACH STATE
% Probabilities are given in P_fullmtx (non-zero for feasible states)
trP_fullmtx=P_fullmtx';
c_state=trP_fullmtx(:); %Get probabilities for all states

c_state(c_state==0)=[]; %Remove zero probability states
%c_state=ones(size(Q,2),1);
  
  %Created LP matrices and vectors.
  
  %----- START ITERATIVE OPTIMIZATION AND BASIS FUNCTION GENERATION -------
  %Initialize counter for number of basis functions
  k=1;
  
  %Set up plotting 
  figure
  plot(cost, '*'); %Plot correct costs
  hold on;
  legend('Actual Cost','Approximate Cost');
  
  %Store approximation error for using k bases (2-NORM)
  approx_err=[];
  %Store approximation
  approx=[];
  
  %Loop to iterate over cost approximation and basis function generation
  while k<=R %Add maximum possible number of bases
      %A) Get kth approximate dual
      %Solve optimization to get vector...
      cvx_begin
        params.OptimalityTol = tolerance; %Set tolerance
        variable r_fit(size(Phi,2))
        dual variables d
        %dual variables d2
        %dual variables d3
        maximize( c_state'*Phi*r_fit ) %<---------------- Phi CHANGES each iteration at end of iteration
        subject to
            %d1: norm(phi_k,1) <= 1e3
            d: Q*Phi*r_fit <= b
            %d3: [Phi phi_k]*[r_fit;1] >= -1e-3
      cvx_end
      optD=d;
      
      %B) Plot; store error
      %Plot updated approximate costs
      plot(Phi*r_fit, '*');
      %Store next approximation error (hopefully LOWER)
      approx_err=[approx_err;norm(cost-Phi*r_fit,2)];
      %Store approximation too
      approx=[approx,Phi*r_fit];
      
      
      %C) Use to get k-th optimal policy
      
      %Get pi
      GetOptPolicyVect;
      
      %D) Get u(x)
      
          %Get policy matrix
          GetCtrlPolicy_DualLP;

          %Create array of state-ctrl tuples, and use to get indices of
          %states associated with each control value
          
          %REMEMBER: indE is an index for an E_Ind in E_Ind_Vect_p
            StateCtrlTupMtx=[];
            indE_u=[];
            for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
                for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
                    for indL=1:size(E_Ind_MtxALL,2)
                        for p=1:P1*P2
                            if(abs(fullPolicyMtx(E_Ind1,E_Ind2,indL,p)-1)<0.01)
                                StateCtrlTupMtx=[StateCtrlTupMtx;E_Ind1,E_Ind2,indL,p]; %Create array
                                %Get state index-control pairs (for indexing within E_Ind_Vect_p)
                                E_Ind_Vect_p=E_Ind_Vect{p};
                                Lmin_offs_p=Lmin_offs{p};
                                indsE=find(E_Ind_Vect_p==((E_Ind1-1)*N2+E_Ind2));
                                indLmin_offs_p=Lmin_offs_p(indsE);
                                indsL_noOffs=indL-indLmin_offs_p(1); %Remove minimum load index offset
                                indE_u=[indE_u;indsE(indsL_noOffs),p];
                            end
                        end
                    end
                end
            end

      
      g_k=[];
      PF_k=[];
      for i=1:length(indE_u) %FOR EACH STATE POSSIBLE, GIVEN THIS POLICY...
          %E) Get k-th vector g
          g_p=g{indE_u(i,2)};           %Get vector of i-th optimal control
          g_k=[g_k;g_p(indE_u(i,1))];   %Add cost of state associated with it

          %F) Get row of matrix PF_k
          P_p=P_mtx{indE_u(i,2)};
          F_p=F{indE_u(i,2)};
          
          %Get row of state associated with said control
          P_p_row=P_p(indE_u(i,1),:);
          %Multiply by F_p to get PF_k row
          PF_k_row=P_p_row*F_p;
          
          %G) Add in to create matrix PF_k
          PF_k=[PF_k;PF_k_row];
      end
      
      %H) Get phi_k+1 (next basis vector)
      phi_k=g_k+DISCOUNT*PF_k*Phi*r_fit;

      %I) Add next basis to Phi
      Phi=[Phi phi_k];

      k=k+1;
  end
  
  title(strcat('Testing Cost Approximation Using', R+1,'-Bases Fit'));
  hold off
  
  figure
  plot(approx_err/norm(cost,2));
  xlabel('Number of Basis Vectors Added');
  ylabel('2-Norm Error in Cost Approximation');
  title('Normalized Cost Approximation Error (Default Test)');
  
  optD = d; %Get vector of FINAL dual
  cost=Phi*r_fit; %Get FINAL approximated cost
  
  %Format FINAL cost vector into E1xE2 matrices (one for each value of load)
  ConvCosts=FormatCostVect(cost);
