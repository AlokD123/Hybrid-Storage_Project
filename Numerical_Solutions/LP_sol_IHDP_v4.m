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
%FullState_Ind_Vec=[];
numAdmissibleLoads=0; %Count number of admissible load values for a given energy state (for UNIFORM DISTRIBUTION)

%F=zeros(M*N1*N2,M*N1*N2);         %Transition matrix, to get next states

%P_fullmtx=zeros(N1*N2,M);         %Full matrix of load probabilities given state
%P=zeros(M*N1*N2,M*N1*N2);         %P matrix, containing probabilties of next states 

%PF=zeros(M*N1*N2,M*N1*N2,P1*P2);  %PF-matrices, product of P and F (one for each value of 1<=p<=P)
%g=zeros(M*N1*N2,P1*P2);           %stage cost (g) vectors (one for each value of 1<=p<=P)
PF={};
P=[];

indL_Feas=[]; %Vector of feasible demands for ONE GIVEN combination of x and u

unrepNextE_Inds=[]; %List of unrepeated nextE_Ind values

%PART A: SET UP MATRICES
%For each possible control...
  for D1=0:MAX_DISCHARGE(1)
    for D2=0:MAX_DISCHARGE(2)
        %Map control to control index
        D1_Ind=D1+1; D2_Ind=D2+1;
        %Get combination #(p)
        p=D2_Ind+P2*(D1_Ind-1);
        
        indCount=0; %Index for feasible state #, for a given value of p
        %Index for row in E-state matrix
        %rowInd_Emtx=p-1; %Reset when restarting from top to add more states for next value of p
        
        %Flag for updating the last stored nextE_Ind value (see below)
        %boolUpdTemp=1; 
        
        %For each state at an iteration...
        for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
            for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
                %Map state index to state
                E1=E_MIN(1)+(E_Ind1-1);
                E2=E_MIN(2)+(E_Ind2-1);
                
                %Get index of current state energies in vector of state energies
                E_Ind=(E_Ind1-1)*N2+E_Ind2;
                
                if(D1>E1 || D2>E2)  %If discharge too high for state, IGNORE
%                     nextE_Ind_Vect=[nextE_Ind_Vect;-1*ones(M,1)]; %No next state index
%                     P_fullmtx(E_Ind,:)=0; %No probable next state
%                     for indL=1:M    %Ignore constraints
%                         g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=INF_COST; %Because probability of transition is zero, have set constraint to ARBITRARY sentinel value
%                     end
                else
                    %Index row in E-state indices mtx (for feasible E-state) same as VALUE of E-state index
                    rowInd_Emtx = E_Ind;
                    %For each perturbation at the CURRENT time
                    for indL=1:(MAX_LOAD-MIN_LOAD+1)
                        %Map index to value of load
                        L=indL+MIN_LOAD-1;

                        %STEP 0
                        %Calculate the state these values of u and w will lead to, even if
                        %impossible...
                        [nextE1,nextE2]=optNextStateLimited(E1,E2,D1,D2,L);

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
                                  E_Ind_MtxALL(rowInd_Emtx,indL)=E_Ind;
                                    
                                  %Map state to state index, to find cost of next state based on its index
                                  nextE_Ind1=round(nextE1-E_MIN(1)+1);
                                  nextE_Ind2=round(nextE2-E_MIN(2)+1); 

                                  %STEP 2
                                  %Get index of next state energy in vector of state energies
                                  nextE_Ind=(nextE_Ind1-1)*N2+nextE_Ind2;
                                  %Add next state energy index to vector of FEASIBLE next state energies
                                  nextE_Ind_Vect_p=[nextE_Ind_Vect_p;nextE_Ind];
                                  
                                  
                                  %STEP
                                  %Create augmented vector containing current E-states - EXCLUDING those nextly infeasible - AND ALSO next E-states
                                  %ONLY for load=0 
                                  if(length(nextE_Ind_Vect_p)==1 && E_Ind==nextE_Ind) %If first state being added, and not differing...
                                      aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %By default, add to augmented vector
                                  else                            %ALL OTHER CASES
                                      if(above_E_Ind==E_Ind) %If repeating E_Ind...
                                          %if(boolUpdTemp) 
                                          unrepNextE_Inds=[unrepNextE_Inds; nextE_Ind]; %Store nextE_Ind value for later (whether currently infeasible OR POSSIBLY NOT)
                                            %boolUpdTemp=0; %Don't update temp until added this state to aug vector
                                          %end
                                          aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %Just add CURRENT E-state to augmented vector
                                      else
                                          %Else, if adding new E_Ind value (not repeating for load)...
                                          %For each next state not (/not yet) in aug_nextE_Ind_Vect_p...
                                          for unrepI=unrepNextE_Inds'
                                              if ~any(nnz(aug_nextE_Ind_Vect_p==unrepI)) %If NOT in state space of current E-states...
                                                aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p; unrepI]; %Insert in between
                                              end
                                          end
                                          aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p;E_Ind]; %Add regular state index at end, as usual
                                          %boolUpdTemp=1; %Resume updating temp
                                          unrepNextE_Inds=nextE_Ind; %Restart list of unrepeated nextE_Ind values
                                      end
                                  end
                                  above_E_Ind=E_Ind; %Update "above" E_ind value to current position
                                  
                                  
                                  
                                  %Add indL to list of FEASIBLE loads for this combination of u and x
                                  indL_Feas=[indL_Feas;indL];
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
                        

                        %STEP 3
                        %Create full probability matrix
                        %SET DISTRIBUTION: UNIFORM <------------------------------------------- ************
                        if(nextE_Ind~=-1) %If this state leads to a feasible next state...
                            %Put sentinel value to indicate as such (for given control)
                            P_fullmtx(E_Ind,indL)=2;
                            %If load admissible, increment count of admissible loads
                            numAdmissibleLoads=numAdmissibleLoads+1;

                            %STEP 4
                            %Create p-th vector g, for constraint
                            %g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=CtrlCost(D1,D2,L);
                            %g(indCount,p)=CtrlCost(D1,D2,L);
                            gVec_p(indCount)=CtrlCost(D1,D2,L); %Cost of stage is given by CtrlCost

                        else %Else if infeasible next state...
                            %DO NOTHING
                            %Set 0 probability (for given control)
                            %P_fullmtx(E_Ind,indL)=0;
                            %Ignore constraint
                            %g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=INF_COST; %Because probability of transition is zero, have set constraint to ARBITRARY sentinel value
                        end
                    end

                    %Fill in row of full matrix with UNIFORM DISTR. values, after finding number of feasible load values
                    for indL_Vect=indL_Feas
                        if(P_fullmtx(E_Ind,indL_Vect)==2)
                           P_fullmtx(E_Ind,indL_Vect)=1/numAdmissibleLoads; %Uniform probability
                        end
                    end
                    %Reset feasible loads count, for subsequent energy state
                    numAdmissibleLoads=0;
                    %Reset list of feasible loads (next state)
                    indL_Feas=[];
                end
            end
        end
    
        %For P_fullmtx entries with no feasible CURRENT state, take NEXT
        %state to be feasible in the no-load case(i.e. P(wk=1|lambda_inf)=1)
        if(~isempty(P_fullmtx(all(P_fullmtx==0,2),:)))
            emptyRows=all(P_fullmtx==0,2);
            for rowNum=find(emptyRows==1)'
                P_fullmtx(rowNum,:)=[1, zeros(1,size(P_fullmtx,2)-1)];
            end
        end
    
    %At end of p-th value cycle, if remaining unadded states, add to end
    for unrepI=unrepNextE_Inds'
        if ~any(nnz(aug_nextE_Ind_Vect_p==unrepI)) %If NOT in state space of current E-states...
            aug_nextE_Ind_Vect_p=[aug_nextE_Ind_Vect_p; unrepI]; %Insert in between
        end
    end
    %Also, exclude from the augmented vector states that are nextly infeasible
    nextlyInfE=~ismember(aug_nextE_Ind_Vect_p,nextE_Ind_Vect_p);
    aug_nextE_Ind_Vect_p(nextlyInfE)=[];
    
    unrepNextE_Inds=[]; %Restart list of unrepeated nextE_Ind values
        
    g{p}=gVec_p';
    E_Ind_Vect{p}=E_Ind_Vect_p;
    nextE_Ind_Vect{p}=nextE_Ind_Vect_p;
    aug_nextE_Ind_Vect{p}=aug_nextE_Ind_Vect_p;
    
    
    
    %STEP 6
    %Create P matrix: select rows corresponding to components in nextE_Ind_Vect
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
    
    %Multiply to get p-th PF matrix
    PF{end+1}=P;

    
    %Reset matrices/vectors
    nextE_Ind_Vect_p=[];
    E_Ind_Vect_p=[];
    aug_nextE_Ind_Vect_p=[];
    
    numAdmissibleLoads=0;
    
    %F=[];           %zeros(M*N1*N2,M*N1*N2);
    P_fullmtx=[];   %zeros(N1*N2,M);
    P=[];           %zeros(M*N1*N2,M*N1*N2);
    
    end
  end
  
  
  %STEP : Construct vector of ALL FEASIBLE energies
  E_Ind_VectALL=[];
  for row=1:size(E_Ind_MtxALL,1)
      nnzRow=nnz(E_Ind_MtxALL(row,:));
      E_Ind_Mtx_nzRow=E_Ind_MtxALL(row,1:nnzRow);
      E_Ind_VectALL=[E_Ind_VectALL; E_Ind_Mtx_nzRow'];
  end
  
  
  %STEP : Construct each F matrix
  for p=1:P1*P2
      aug_nextE_Ind_Vect_p=aug_nextE_Ind_Vect{p};
      %Index COLUMN of F matrix by ROW number of E_Ind_VectALL
      row=1; %Reset row being checked in E_Ind_VectALL to start when start on next E_Ind vector
      
      %Go through next E-state index vector for current value of p...
      for r=1:length(aug_nextE_Ind_Vect_p)
          %If next state is currently infeasible...
          if aug_nextE_Ind_Vect_p(r)<E_Ind_VectALL(row) %(i.e. NOT continuously increasing in augmented vector)
             row=1; %Restart from beginning of E_Ind_VectALL to find the state <---------------------------------- ASSUMING ONLY 1 new currently infeasible state!
          end
          while(E_Ind_VectALL(row)~=aug_nextE_Ind_Vect_p(r)) %While not reached mapping column in F (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          F_p(r,row)=1; %Once reached, map
          row=row+1;  %Start at next column in F next time <----------------------------------------------------- ASSUMING continuously increasing in augmented vector!!!
      end
      %Add extra zeros at end to ensure dimensions of F_p and E_Ind_VectALL match
      F_p(:,row:length(E_Ind_VectALL))=0;
      
      F{p}=F_p;
      F_p=[]; %Reset
      
      PF{p} = PF{p}*F{p}; %Finish PF matrices
  end
  
  %STEP : Construct each G matrix
  for p=1:P1*P2
      E_Ind_Vect_p=E_Ind_Vect{p};
      %Index COLUMN of G matrix by ROW number of E_Ind_VectALL
      row=1; %Reset row being checked in E_Ind_VectALL to start when start on next E_Ind vector
      
      %Go through E-state index vector for current value of p...
      for r=1:length(E_Ind_Vect_p)
          while(E_Ind_VectALL(row)~=E_Ind_Vect_p(r)) %While not reached mapping column in G (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          G_p(r,row)=1; %Once reached, map
          row=row+1;  %Start at next column in G next time (continuously increasing)
      end
      %Add extra zeros at end to ensure dimensions of G_p and E_Ind_VectALL match
      G_p(:,row:length(E_Ind_VectALL))=0;
      
      G{p}=G_p;
      G_p=[]; %Reset
  end
  
  c;
  %PART B: OPTIMIZATION
  %Created LP matrices and vectors.
  %Run optimization problem, and find primal as well as dual.
  cvx_begin
    grbControl.LPMETHOD = 1; % Use dual simplex method
    variable cost(length(E_Ind_VectALL))
    dual variables d{P1*P2}
    maximize( sum(cost) )
    subject to
        for p=1:P1*P2
            %d{p} : (eye(length(PF))-DISCOUNT*PF(:,:,p))*cost <= g(:,p)
            d{p} : (G{p}-DISCOUNT*PF{p})*cost <= g{p} %<-------- HOW TO ENSURE PF{p} is SQUARE for ALL p???
        end
  cvx_end
  %Get vector of optimal dual from cell array form
  optD = cell2mat(d);
  
  %Format cost vector into E1xE2 matrices (M matrices, for each value of load)
    for(i1=0:M*N2:M*N2*(N1-1))
        for(j=0:M-1)
          for(ind=(1+i1+j):M:(M*(N2-1)+i1+(j+1)))
              if(mod(ind,M*N2)==0)
                ind2=(M*N2-1-j)/M+1;
              else
                ind2=(mod(ind,M*N2)-1-j)/M+1;
              end
              ConvCosts(i1/(M*N2)+1,ind2,j+1)=cost(ind);
          end
        end
    end
    
    %PART C: STATIONARY POLICY
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