%LP solution of IHDP (Value Iteration) for Hybrid Storage optimization
%warning('off', 'Octave:possible-matlab-short-circuit-operator');
clearvars

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
DISCOUNT=0.99;


%Definitions
M=MAX_LOAD-MIN_LOAD+1;
N1=(E_MAX(1)-E_MIN(1)+1);
N2=(E_MAX(2)-E_MIN(2)+1);
P1=MAX_DISCHARGE(1)+1;
P2=MAX_DISCHARGE(2)+1;

%Initialization
%E_Ind_Vec=[];%FullState_Ind_Vec=[];
nextE_Ind_Vect=[]; %Vector of next state energies
numAdmissibleLoads=0; %Count number of admissible load values for a given energy state (for UNIFORM DISTRIBUTION)

F=zeros(M^2*N1*N2,M^2*N1*N2);       %Transition matrix, to get next states

P_fullmtx=zeros(N1*N2,M);         %Full matrix of load probabilities given state
P=zeros(M*N1*N2,M^2*N1*N2);       %P matrix, containing probabilties of next states 

%Create "identity" matrix, duplicating elements of state energy vector, for
%mapping current energy state to combination of new state and next load
v=ones(M,1); %Vector of ones along diagonal, which is repeated to repeat values of state
Id = kron(eye(M*N1*N2),v);

PFId=zeros(M*N1*N2,M*N1*N2,P1*P2);  %PFId-matrices, product of P, F, and Id (one for each value of 1<=p<=P)
g=zeros(M*N1*N2,P1*P2);           %stage cost (g) vectors (one for each value of 1<=p<=P)

%For each possible control...
  for D1=0:MAX_DISCHARGE(1)
    for D2=0:MAX_DISCHARGE(2)
        %Map control to control index
        D1_Ind=D1+1; D2_Ind=D2+1;
        %Get combination #(p)
        p=D2_Ind+P2*(D1_Ind-1);
        
        %For each state at an iteration...
        for E_Ind1=1:(E_MAX(1)-E_MIN(1)+1)
            for E_Ind2=1:(E_MAX(2)-E_MIN(2)+1)
                %Map state index to state
                E1=E_MIN(1)+(E_Ind1-1);
                E2=E_MIN(2)+(E_Ind2-1);
                
                if(D1>E1 || D2>E2)  %If discharge too high for state, no constraints
                    nextE_Ind_Vect=[nextE_Ind_Vect;-1*ones(M,1)]; %No next state index
                    P_fullmtx(:,:)=0; %No probable next state
%                     for indL=1:M    %Ignore constraints
%                         g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=inf; 
%                     end
                else
                    %Get index of current state energies in vector of state energies
                    E_Ind=(E_Ind1-1)*N2+E_Ind2;

                    %For each perturbation at the CURRENT time
                    for indL=1:(MAX_LOAD-MIN_LOAD+1)
                        %Map index to value of load
                        L=indL+MIN_LOAD-1;

                        %STEP 1
                        %Calculate the state these values of u and w will lead to, even if
                        %impossible...
                        [nextE1,nextE2]=optNextStateLimited(E1,E2,D1,D2,L);

                        %If next state is amongst those achievable with a given perturbance....
                        if(nextE1<=E_MAX(1) && nextE1>=E_MIN(1))
                            if(nextE2<=E_MAX(2) && nextE2>=E_MIN(2))
                                %IF meeting following conditions: (C_MAX and C_MIN)
                                %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
                                if(~((D1+D2-L)<0||(D1+D2-L)>MAX_CHARGE(2)))
                                  %Map state to state index, to find cost of next state based on its index
                                  nextE_Ind1=round(nextE1-E_MIN(1)+1);
                                  nextE_Ind2=round(nextE2-E_MIN(2)+1); 

                                  %%NEW
                                  %STEP 2
                                  %Get index of next state energies in vector of state energies
                                  nextE_Ind=(nextE_Ind1-1)*N2+nextE_Ind2;
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

                        %Add next state energy index to vector of feasible next state energies
                        nextE_Ind_Vect=[nextE_Ind_Vect;nextE_Ind];

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
                            %g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=CtrlCost(D1,D2,L); %Cost of stage is given by CtrlCost

                        else %Else if infeasible next state...
                            %Set 0 probability (for given control)
                            P_fullmtx(E_Ind,indL)=0;
                            %g(indL+M*(E_Ind2-1+N2*(E_Ind1-1)),p)=inf; %Ignore constraint
                        end
                    end

                    %Fill in row of full matrix with UNIFORM DISTR. values, after finding number of feasible load values
                    for indL=1:(MAX_LOAD-MIN_LOAD+1)
                        if(P_fullmtx(E_Ind,indL)==2)
                           P_fullmtx(E_Ind,indL)=1/numAdmissibleLoads; %Uniform probability
                        end
                    end
                    %Reset feasible loads count, for subsequent energy state
                    numAdmissibleLoads=0;
                end
            end
        end
    
    %STEP 5
    %Create duplicated next state index vector, to determine next state for
    %each random perturbation, @ each SINGLE current energy state
    dup_nextE_Ind_Vect=[];
    for i=1:length(nextE_Ind_Vect)
        dup_nextE_Ind_Vect=[dup_nextE_Ind_Vect;repmat(nextE_Ind_Vect(i),M,1)];%Repeat each new state index M times
    end
    
    %Construct F matrix
    for r=1:M^2*N1*N2
      boolRowFull=0; %Flag to achieve 1 element per row of F. Set flag off initially.
      indDupVect=r; %Get index for next element in dup_nextE_Ind_Vect to check as r-th row, since next state index is constant for all elements in row
      
      for c=1:M^2*N1*N2
        if(dup_nextE_Ind_Vect(indDupVect)==-1) %If find infeasible next state...
            F(r,c)=0; %Do not consider state
        else  %Else, fill in row with 1's so as to map state to next state
            if (mod(r,M)==0) %Special case for mod function w/ 1-indexing
                Cond=(c==(dup_nextE_Ind_Vect(indDupVect))*M);                          %<----------------------- TO CHECK!!
            else
                Cond=(c==((dup_nextE_Ind_Vect(indDupVect)-1)*M+mod(r,M)));             %<----------------------- TO CHECK!!
            end
            %^^ General next state mapping condition. 1 if mapping.
            if (Cond &&(~boolRowFull)) %Full condition. If generally mapping and not full, set mapping (fill 1).
                F(r,c)=1; 
            end
        end
        
      end
    end

    %STEP 6
    %Create P matrix: take select rows corresponding to components in nextE_Ind_Vect
    for r=1:M*N1*N2
        c=(r-1)*M+1; %Get column number of next row of probabilities (given state)
        Ind_nextE=nextE_Ind_Vect(r);    %Get index of state stored in r-th row of nextE_Ind_Vect (i.e. the next energy state)
        
        if(Ind_nextE==-1) %If find infeasible next state...
            P(r,c:(c+M-1))=zeros(1,M); %Fill in row r with 0's
        else
            P(r,c:(c+M-1))=P_fullmtx(Ind_nextE,:); %Else, fill in row r with probabilities associated with Ind_nextE-th row of P_fullmtx
        end
    end
    
    %Multiply to get p-th PFId matrix
    PFId(:,:,p)=P*F*Id;

    
    %Reset matrices/vectors
    nextE_Ind_Vect=[];
    numAdmissibleLoads=0;
    F=zeros(M^2*N1*N2,M*N1*N2);
    P_fullmtx=zeros(N1*N2,M);
    P=zeros(M*N1*N2,M^2*N1*N2);
    dup_nextE_Ind_Vect=[];
    
    end
  end
  
  %Optimization
  
  
%   %%%%% EXTRA... Not needed?
%   
% %Add to vector of FULL states, (E1,E2,L).... Will map
% %1-1 to vector of costs
% FullState_Ind_Vec=[FullState_Ind_Vec;E_Ind1,E_Ind2,indL];
%   
%   %Find probability for each possible value of perturbation (w)
%   %NOTE: this is a perturbation at the NEXT time (leading to an expected cost-to-go at the CURRENT time)
%   for indNextL=1:(MAX_LOAD-MIN_LOAD+1)
%       %Map index to value of load
%       NextL=indNextL+MIN_LOAD-1;
%       P(indNextL,indL,)              %% NOTE: using wj instead of f(xi,u,wj) in conditional probability because not need for next state indexing (load indexing is equivalent)
%       %Pi=[Pi ];
%   end
% 
%   
% %Add to vector of energy states, (E1,E2)
% E_Ind_Vec=[E_Ind_Vec;E_Ind1,E_Ind2];