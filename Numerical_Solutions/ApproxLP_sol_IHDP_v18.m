%ALP solution of IHDP (Value Iteration) for Hybrid Storage optimization
% ** ADDED REGENERATIVE BRAKING **
% COMBINED CONTROLS (default mutual exclusion)
% **** NEW STAGE COST ****
%Realistic sizing

clearvars -except cost approx_err E_MAX max_E_SIZE min_E_SIZE minCost max_E1 max_E2 optE_SIZE c1 c2 vectS_netOptVal PF_opt_mtx g_opt_vect P_opt_mtx size_iter PF_opt g_opt Exp_CostToGo feasStatesArr_size optVal_size optCost_size g_opt_mtx Exp_CostToGo_mtx totCost INFCOST size1_mult size2_mult RES_E1 RES_E2 RES_L RES_U1 SCALE_BATT SCALE_SC SCALE_L Phi_size maxE_stepSize E1_counter E2_counter maxE_stepSize_E1 maxE_stepSize_E2 mult_cost_idx cost_mult_1 cost_mult_2 optRatio optCosts c1_fact c2_fact optSize opt_E_size P0;

global E_MIN;
E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
global E_MAX; 
%{
E_MAX=[50000,250]; %Maximum energy to be stored (upper bound)        %Comment out when running storage sizing
%E_MAX=[10,5];
%}

%Solver tolerance
tolerance=1e-6;

%% Input: initial state, horizon, number of bases in approximation
%Initial stored energy (user-defined)
%Must be between MIN_STATE and MAX_STATE
E1_INIT=E_MAX(1); 
E2_INIT=E_MAX(2);

R=RES_E2*E_MAX(2);   %MAXIMUM order of extra polynomial bases added by iteration (TOTAL MUST BE LESS THAN NUMBER OF FEASIBLE STATES)
MAX_STEPS=10; %MAXIMUM number of groups in state aggregation

%% Model setup
global MAX_CHARGE; global MAX_DISCHARGE;
mult=1;

%{
MAX_CHARGE=[1*2*size1_mult;15*2*size2_mult]; %Maximum charging of the 1) battery and 2) supercap
MAX_DISCHARGE=[1*2*size1_mult;15*2*size2_mult]; %Maximum discharging of the 1) battery and 2) supercap
%}
%
MAX_CHARGE=mult*[1,15];                                                           %Comment out when running storage sizing
MAX_DISCHARGE=mult*[1,15];
%}

global MIN_LOAD; global MAX_LOAD; global SCALE_L;
%
MIN_LOAD=-(MAX_CHARGE(1)+MAX_CHARGE(2)); %Maximum regenerative energy expected
MAX_LOAD=MAX_DISCHARGE(1)+MAX_DISCHARGE(2); %Maximum load expected
%}
SCALE_L=(4.5+7.5)/MAX_LOAD; %Wh/gridpt

%{
MIN_LOAD=-mult*2;
MAX_LOAD=mult*2;
%}

MAX_NUM_ZEROS=3; %Maximum number of zero load counts before end sim

global ALPHA_C; global ALPHA_D; global BETA; global K; global nu;
ALPHA_C=[0.87;0.84]; %Efficiency of charging
ALPHA_D=[0.87;0.84]; %Efficiency of discharging
BETA=[1;1];       %Storage efficiency
K=[1e3;1e3];      %Weighting factors for D1^2 and C1^2 costs
nu=1;  %Cost for high state of charge
PERFECT_EFF=0;

%Discounted infinite horizon problem
global DISCOUNT; %Discount factor
DISCOUNT=0.99;


%% Definitions
%Resolution of battery charge/dischage
global RES_U1;
%Resolution of state variables (GRID RESOLUTION)
global RES_E1;
global RES_E2;
global RES_L;
%Resolution of demand DURING SIMULATION ("continuous" demands)
global resL_Mult;
resL_Mult=RES_L;
%{
RES_U1=1;
RES_E1=1/1000;
RES_E2=1/5;
RES_L=1;
resL_Mult=10;  %NOT necessarily equal to grid resolution for DP             %<----CHANGE THIS TO CHANGE THE SIMULATION RUNTIME!!!!
%}

global N2; global N1; global M; global P0;

M=round2even(RES_L*(MAX_LOAD-MIN_LOAD))+1;
N1=RES_E1*(E_MAX(1)-E_MIN(1))+1;
N2=RES_E2*(E_MAX(2)-E_MIN(2))+1;

P0=RES_U1*(MAX_CHARGE(1)+MAX_DISCHARGE(1))+1;

global epsilon; global epsilon2; global epsilon3;
global epsilon4; global epsilon5; global gamma; global INF_Q
epsilon=0.01; %Next state off grid rounding tolerance
epsilon2=0.0001; %Off grid state comparison tolerance
epsilon3=0.01; %On-grid optimal control search tolerance
epsilon4=1e-5; %On-grid optimal q search tolerance
epsilon5=0.01; %Number of repeated q-values count tolerance

gamma=1e-2; %Regularization term weighting factor
INF_Q=Inf; %Define infeasible Q-value to be infinite (so cannot choose as minimal)


%% Initialization
E_Ind_Vect_p=[];      %Vector of current state energies
E_Ind_L_Vect_p=[];    %Array of vectors of current state energies and associated feasible demands FOR EACH energy state
nextE_Ind_Vect_p=[];  %Vector of next state energies
nextE_Ind_L_Vect_p=[]; %Array of vectors of next state energies and associated feasible demands FOR EACH energy state
aug_nextE_Ind_Vect_p=[]; %Augmented vector containing current state energies and next energies for currently infeasible states
aug_Vect_Ls_p=[]; %Same, but also including associated load values in each state

offGrdNxtE1E2_p=[]; %Array mapping single index to linear next state index, for states OFF THE GRID
numLoads_OffGrd_p=[]; %Vector of number of admissible loads in next states that are OFF THE GRID

numL_OffGrd=0; %Count number of admissible load values for a given NEXT energy state

global E_Ind_MtxALL; %Matrix of all states (0 if infeasible)
global E_Ind_VectALL; %Vector of all feasible states
global E_VectALL_Ls; %Vector of all associated loads
global feasInf_E_VectALL_Ls; %Vector of both feasible associated loads and infeasible ones
global CostMtx; %Matrix with 3 columns: a) E-state, b) load, c) associated cost. Used for approximation
CostMtx=[];
E_Ind_MtxALL=[];
E_Ind_VectALL=[];
E_VectALL_Ls=[];
feasInf_E_VectALL_Ls=[];

P_mtx={};   %Array of P matrices
P=[];       %Current P matrix
PF={};      %Array of P*F matrices

F_p=[]; G_p=[]; %F and G matrices, for mapping states

fbleEStates_p=[]; %Vector of indicating which E-statea are feasible for given net discharge (value of p)
Lmin_offs_p=[]; %Vector of minimum load offsets for each E-state, to create CORRECT MAPPING in G matrix
E_Ind_Mtx_p=[]; %Matrix of E_Ind_MtxALL values, but for EACH value of p
P_fullmtx=[];   %Matrix of all probabilities

indL_Feas=[]; %Vector of feasible demands for ONE GIVEN combination of x and u
feasStates=[]; %List of all feasible states (E1,E2,L), no repeats
feasE1s_p=[]; feasE2s_p=[]; feasLs_p=[]; %Vectors of feasible states, for each p
feasStates_p=[]; %Array of adjoint vectors, for each p

p_max=0;    %Maximum number of controls to consider, initialized at 0


Phi=[]; %Design matrix, for cost approximation
margApproxErr=zeros(N1,1); %Error in each state marginalized, initialized at 0s
c_state=[];     %Vector of state-relevance weightings

gVec_p=[];
%% PART A: SET UP MATRICES
%For each possible control...
  for U1_Ind=1:P0
      %Map control to control index
        U1=(U1_Ind-1)/RES_U1-MAX_CHARGE(1);
        %Get combination #(p)
        p=U1_Ind;

        indCount=0; %Index for feasible state #, for a given value of p

        %For each state at an iteration...
        for E_Ind1=1:N1
            for E_Ind2=1:N2
                %Map state index to state
                E1=E_MIN(1)+(E_Ind1-1)/RES_E1;
                E2=E_MIN(2)+(E_Ind2-1)/RES_E2;

                %Get index of current state energies in vector of state energies
                E_Ind=(E_Ind1-1)*N2+E_Ind2;

                boolFbleState=0; %Flag for checking if state is feasible. Flag=1 only if at least one demand leads to a feasible state&next-state

                %If net discharge is too low (negative) or high for battery
                if(U1>(BETA(1)*E1-E_MIN(1))*ALPHA_D(1) || U1<(BETA(1)*E1-E_MAX(1))/ALPHA_C(1))
                    %IGNORE, but still account for state being counted by inserting row of zeros
                    rowInd_Emtx = E_Ind;
                    E_Ind_Mtx_p(rowInd_Emtx,:)=0;
                else
                    %Index row in E-state indices mtx (for feasible E-state) same as VALUE of E-state index
                    rowInd_Emtx = E_Ind;
                    %Store number of feasible E-states in which can  apply control
                    fbleEStates_p=[fbleEStates_p;E_Ind1,E_Ind2];

                    %For each perturbation at the CURRENT time...
                    for indL=1:M
                        %Map index to value of load
                        L=MIN_LOAD+(indL-1)/RES_L;

                        %STEP 0
                        %Calculate the state these values of u and w will lead to, even if
                        %impossible...
                        [nextE1,nextE2]=optNextStateLimited_v3(E1,E2,U1,L);

                        %If next state is amongst those achievable with a given perturbance....
                        if(nextE1<=E_MAX(1) && nextE1>=E_MIN(1))
                            if(nextE2<=E_MAX(2) && nextE2>=E_MIN(2))
                                %IF meeting following conditions: (C_MIN and C_MAX)
                                %1) net supply (discharging) never below demand, 2) not charging cap. too quickly
                                if(~( (-MAX_CHARGE(2)+U1)>L || L>(MAX_DISCHARGE(2)+U1) ))

                                    boolFbleState=1; %State&next-state are feasible for this demand L

                                  %Count the number of feasible states for a given set of controls (D1,D2)
                                  indCount=indCount+1; %... and use as an index

                                  %STEP 1: create vector and matrix of FEASIBLE state energies for each load
                                  %Add state energy index to vector for current value of p (D1,D2 combo)
                                  E_Ind_Vect_p=[E_Ind_Vect_p;E_Ind];
                                  E_Ind_L_Vect_p=[E_Ind_L_Vect_p;E_Ind,L];
                                  %Add state energy index to matrix of ALL FEASIBLE energies
                                  %DO NOT RESET at end. Will overwrite with same values (and add) each time, which is ok.
                                  E_Ind_MtxALL(rowInd_Emtx,indL)=E_Ind;
                                  E_Ind_Mtx_p(rowInd_Emtx,indL)=E_Ind;

                                  %Map state to state index, to find cost of next state based on its index
                                  if(abs((RES_E1*(nextE1-E_MIN(1))+1)-round(RES_E1*(nextE1-E_MIN(1))+1)<epsilon))
                                      nextE_Ind1=round(RES_E1*(nextE1-E_MIN(1))+1);  %IF within error bound, round
                                  else
                                      nextE_Ind1=RES_E1*(nextE1-E_MIN(1))+1;         %.. Otherwise, will interpolate
                                  end
                                  if(abs((RES_E2*(nextE2-E_MIN(2))+1)-round(RES_E2*(nextE2-E_MIN(2))+1))<epsilon)
                                      nextE_Ind2=round(RES_E2*(nextE2-E_MIN(2))+1);
                                  else
                                      nextE_Ind2=RES_E2*(nextE2-E_MIN(2))+1;
                                  end

                                  %STEP 2: create vector of next state energies for each load
                                  %Get index of next state energy in vector of state energies
                                  nextE_Ind=(nextE_Ind1-1)*N2+nextE_Ind2;
                                  %Add next state energy index to vector of FEASIBLE next state energies
                                  nextE_Ind_Vect_p=[nextE_Ind_Vect_p;nextE_Ind];


                                  %STEP A: create array of next state indices that are off-grid
                                  %Mapping from 2 indices to linear index
                                  if round(nextE_Ind)~=nextE_Ind %If NOT INTEGER...... i.e. EITHER component of next state is off grid...
                                      offGrdNxtE1E2_p=[offGrdNxtE1E2_p;nextE_Ind1,nextE_Ind2,nextE_Ind];
                                  end


                                  %STEP 3: determine feasible loads
                                  %Add indL to list of FEASIBLE loads for this combination of u and x
                                  indL_Feas=[indL_Feas;indL];

                                  %STEP 4: Create boolean array of all FEASIBLE states
                                  feasStates(E_Ind1,E_Ind2,indL)=1;
                                  %Also create separate lists just for current control
                                  feasE1s_p=[feasE1s_p;E1]; feasE2s_p=[feasE2s_p;E2]; feasLs_p=[feasLs_p;L];
                                  feasStatesArr_p=[RES_E1*(feasE1s_p-E_MIN(1))+1,RES_E2*(feasE2s_p-E_MIN(2))+1,RES_L*(feasLs_p-MIN_LOAD)+1];%Adjoin into array

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
                            gVec_p(indCount)=CtrlCost_Modified_v3(U1,L,E1); %Cost of stage is given by CtrlCost
                        else %Else if infeasible next state...
                            %DO NOTHING
                        end
                    end

                    %If NO LOAD LED TO FEASIBLE STATE, set that row in E_Ind_Mtx_p to zeros
                    if(~boolFbleState)
                        E_Ind_Mtx_p(rowInd_Emtx,:)=0;
                    end

                    %Reset list of feasible loads (next state)
                    indL_Feas=[];
                end

            end
        end

         %Count number of feasible loads for NEXT E-states off grid (NOT FOR CURRENT ONES)
         prev_Delta2_next=0; %Initialization

      %STEPS B-D
      %Count for ALL possible controls in NEXT STATE...
      for i=1:size(nextE_Ind_Vect_p,1) %For each next state...
          %Maximum load previously achieved by all ctrls (starting point, to not double-count)
          maxL_prev=0; %Reset maximum previous load value

          if round(nextE_Ind_Vect_p(i))==nextE_Ind_Vect_p(i)
              %If next state is on grid, calculate components regularly
              nextE2_Ind=remainder(nextE_Ind_Vect_p(i),N2);
              nextE1_Ind=(nextE_Ind_Vect_p(i)-nextE2_Ind)/N2+1;
              nextE1=(nextE1_Ind-1)/RES_E1+E_MIN(1); nextE2=(nextE2_Ind-1)/RES_E2+E_MIN(2); 
          else
              %Otherwise, if off grid...
              %Find in table of off-grid states %<----------------------------------- ASSUMING ALL OFF-GRID INDICES ARE UNIQUE
              ind_offGrd=find(abs(offGrdNxtE1E2_p(:,3)-nextE_Ind_Vect_p(i))<epsilon2,1);
              if isempty(ind_offGrd)
                  disp("Error, next state off grid not found!");
              end
              %Then get component values
              nextE1=(offGrdNxtE1E2_p(ind_offGrd,1)-1)/RES_E1+E_MIN(1); nextE2=(offGrdNxtE1E2_p(ind_offGrd,2)-1)/RES_E2+E_MIN(2);
          end


          for U1_Ind_next=1:P0
                U1_next=(U1_Ind_next-1)/RES_U1-MAX_CHARGE(1);
                %Get number of FEASIBLE next loads
                %Check excess discharge condition
                if( ~(U1_next>nextE1 || U1_next<(nextE1-E_MAX(1))) )
                    %For each perturbation at the NEXT time...
                    for indL=1:M
                        %Map index to value of load
                        L=MIN_LOAD+(indL-1)/RES_L;                          %Next L is defined to be ON GRID!! (NOT OFF GRID due to control resolution)
                        [next_nextE1,next_nextE2]=optNextStateLimited_v3(nextE1,nextE2,U1_next,L);
                        %Check other conditions
                        if(next_nextE1<=E_MAX(1) && next_nextE1>=E_MIN(1))
                            if(next_nextE2<=E_MAX(2) && next_nextE2>=E_MIN(2))
                                if(~( (-MAX_CHARGE(2)+U1_next)>L || L>(MAX_DISCHARGE(2)+U1_next) ))
                                    %If feasible...
                                    %Get each load in this NEXT energy state
                                    nextE_Ind_L_Vect_p=[nextE_Ind_L_Vect_p;nextE_Ind_Vect_p(i),L];
                                end
                            end
                        end
                    end
                end
          end

          %STORE ONLY UNIQUE LOADS PER NEXT STATE
          nextE_Ind_L_Vect_p=unique(nextE_Ind_L_Vect_p,'rows','stable');

         if round(nextE_Ind_Vect_p(i))~=nextE_Ind_Vect_p(i) %ONLY if off grid..
             %Store count of unique loads for this next state
             numLoads_OffGrd_p=[numLoads_OffGrd_p;nnz(abs(nextE_Ind_L_Vect_p(:,1)-nextE_Ind_Vect_p(i))<epsilon2)];
             %Reset feasible loads count, for subsequent NEXT energy state
             numL_OffGrd=0;
         end
      end    

        %If at least one feasible state, consider this control
        if (~isempty(gVec_p))
            %Store vector data in cell array
            g{p}=gVec_p';
            E_Ind_Vect{p}=E_Ind_Vect_p;
            nextE_Ind_Vect{p}=nextE_Ind_Vect_p;
            fbleEStates{p}=fbleEStates_p;
            Lmin_offs{p}=Lmin_offs_p;
            E_Ind_Mtx{p}=E_Ind_Mtx_p;
            offGrdNxtE1E2{p}=offGrdNxtE1E2_p;
            numLoads_OffGrd{p}=numLoads_OffGrd_p;
            feasStatesArr_ctrl{p}=feasStatesArr_p;
            nextE_Ind_L_Vect{p}=nextE_Ind_L_Vect_p;
            E_Ind_L_Vect{p}=E_Ind_L_Vect_p;

            %Continue testing next control
            p_max=p_max+1;
        else
            %Else, IGNORE

            %Store vector data in cell array
            g{p}=[];
            E_Ind_Vect{p}=[];
            nextE_Ind_Vect{p}=[];
            fbleEStates{p}=[];
            Lmin_offs{p}=[];
            E_Ind_Mtx{p}=[];
            offGrdNxtE1E2{p}=[];
            numLoads_OffGrd{p}=[];
            %feasE1s_ctrl{p}=feasE1s_p; feasE2s_ctrl{p}=feasE2s_p; feasLs_ctrl{p}=feasLs_p;
            feasStatesArr_ctrl{p}=[];
            nextE_Ind_L_Vect{p}=[];
            E_Ind_L_Vect{p}=[];

            %Continue testing next control
            p_max=p_max+1;

        end

        %Reset matrices/vectors
        nextE_Ind_Vect_p=[];
        E_Ind_Vect_p=[];
        gVec_p=[];
        fbleEStates_p=[];
        Lmin_offs_p=[];
        E_Ind_Mtx_p=[];
        offGrdNxtE1E2_p=[];
        numLoads_OffGrd_p=[];
        feasE1s_p=[]; feasE2s_p=[]; feasLs_p=[];
        feasStatesArr_p=[];
        nextE_Ind_L_Vect_p=[];
        E_Ind_L_Vect_p=[];

  end

  %% REMOVE LOADS NOT FEASIBLE FOR ANY E-STATE <----------------------------------------------------------------------------------------------------------------------------------------------
  %{
  minLoad_inc=0;
  i=1;
  while i<=size(feasStates,3) %Iterate through loads
    if all(all(feasStates(:,:,i)==0))
        feasStates(:,:,i)=[]; %Remove
        minLoad_inc=minLoad_inc+1; %Increment minimum load if removed
    end
    i=i+1;
  end
  
  %}

  %%

  %STEP 6: Get MINIMUM load for every E-state in E_VectALL (i.e. baseline for
  %indexing loads for that E-state). NOT SAME AS L_min_offs, which is the
  %minimum offset given a PARTICULAR set of controls.
  
  %MIN_LOAD=MIN_LOAD+minLoad_inc*1/RES_L;
  %M=round(RES_L*(MAX_LOAD-MIN_LOAD)+1);
 
  for i=1:size(feasStates,1)
    for j=1:size(feasStates,2)
        k=1;
        while feasStates(i,j,k)~=1 %Feasible state matrix contains 0s until minimum feasible L
            k=k+1; %So increment till 1
        end
       MINLoad_E_state(j+N2*(i-1),1)=(k-1)/RES_L+MIN_LOAD; %k is now index of minimum L <----------------------------------------------------------------------------------------------------------------------ASSUMING k=1 corresponds to MIN_LOAD!!!!!!!!!!!!!!!
    end
  end

  %STEP 7: Construct vector of ALL FEASIBLE energies, for all controls
  E_Ind_VectALL=[];
  E_VectALL_Ls=[]; %Vector with associated loads
  feasInf_E_VectALL_Ls=[]; %Same, but also including INF where loads are infeasible
  for row=1:size(E_Ind_MtxALL,1)
      nnzRow=nnz(E_Ind_MtxALL(row,:));
      %E_Ind_Mtx_nzRow=E_Ind_MtxALL(row,1:nnzRow);
      E_Ind_VectALL=[E_Ind_VectALL; nonzeros(E_Ind_MtxALL(row,:))];
      %Create vector of all loads per E-state (for all controls)
      Loads_E_state=(find(E_Ind_MtxALL(row,:)>0)-1)/RES_L+MIN_LOAD;
      E_VectALL_Ls=[E_VectALL_Ls;Loads_E_state']; %Append to E_VectALL_Ls
      %numLoads_lt_MINLoad=nnz(~E_Ind_MtxALL(row,1:RES_L*(MINLoad_E_state(row)-MIN_LOAD)+1));
      %numLoads_gt_MAXLoad=nnz(~E_Ind_MtxALL(row,(MINLoad_E_state(row)+nnzRow):size(E_Ind_MtxALL,2)));
      %feasInf_E_VectALL_Ls=[feasInf_E_VectALL_Ls;Inf*ones(numLoads_lt_MINLoad,1);(MINLoad_E_state(row):(MINLoad_E_state(row)+nnzRow-1))';Inf*ones(numLoads_gt_MAXLoad,1);];
  end


  %STEP 8: Create full probability matrix
  %DISTRIBUTION: UNIFORM (unused)
  %(NOTE: actual probabilities are UNUSED!!!)
  for r=1:size(E_Ind_MtxALL,1)
      P_fullmtx(r,:)=E_Ind_MtxALL(r,:)/sum(E_Ind_MtxALL(r,:)); %<----------------For UNIFORM probability, just NORMALIZE rows of feasible states!!
  end

  for p=1:p_max
    if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...

    E_Ind_Vect_p=E_Ind_Vect{p};
    nextE_Ind_Vect_p=nextE_Ind_Vect{p};
    Lmin_offs_p=Lmin_offs{p};
    numLoads_OffGrd_p=numLoads_OffGrd{p};
    nextE_Ind_L_Vect_p=nextE_Ind_L_Vect{p};
    feasStatesArr_p=feasStatesArr_ctrl{p};

    %STEP 9: Create augmented vector containing current E-states - EXCLUDING those nextly infeasible - AND ALSO next E-states
    %(Note: doing after E_Ind_VectALL complete)
    augVectRow=1; %Index row in new augmented vector
    r=1; offGrdNxtE_Idx=1; %Start from beginning
    while r<(length(nextE_Ind_Vect_p)+1) %For each next E-state WITH CURRENT CONTROL COMBO (p)
        if (r~=1) %...IN MOST CASES
            %If next E-state already counted once, do not double-count...
            while r<(length(nextE_Ind_Vect_p)+1) && nnz(abs(nextE_Ind_Vect_p(1:r-1)-nextE_Ind_Vect_p(r))<epsilon2)
                if round(nextE_Ind_Vect_p(r))~=nextE_Ind_Vect_p(r) %If off grid, also skip to subsequent unrepeated number of loads (associated)
                    offGrdNxtE_Idx=offGrdNxtE_Idx+1;
                end
                r=r+1; %Skip to next unrepeated E-state
            end
        end
        if(r~=(length(nextE_Ind_Vect_p)+1))
            %Determine TOTAL number of possible loads for that E-state, given ANY POSSIBLE control used

            %If off grid... STEP E
            if round(nextE_Ind_Vect_p(r))~=nextE_Ind_Vect_p(r)
                numRepNextE=numLoads_OffGrd_p(offGrdNxtE_Idx);
                offGrdNxtE_Idx=offGrdNxtE_Idx+1;
            else %Otherwise...        
                numRepNextE=nnz(E_Ind_VectALL==nextE_Ind_Vect_p(r)); %Number of possible loads is number of times repeated in E_Ind_VectALL
            end

            %Add given E-state to augmented vector that many times (for each load)
            aug_nextE_Ind_Vect_p(augVectRow:(augVectRow+numRepNextE-1),1)=nextE_Ind_Vect_p(r);
            indNextL=find(abs(nextE_Ind_L_Vect_p-nextE_Ind_Vect_p(r))<epsilon2,1); %Get index of first nextE_Ind element in vector with associated loads
            aug_Vect_Ls_p(augVectRow:(augVectRow+numRepNextE-1),1)=nextE_Ind_L_Vect_p(indNextL:(indNextL+numRepNextE-1),2);
            augVectRow=augVectRow+numRepNextE; %Start adding at end next time 
        end
        r=r+1; %Manually increment index in while loop
    end

    %Also, exclude from the augmented vector states that are nextly infeasible
    nextlyInfE=~ismember(aug_nextE_Ind_Vect_p,nextE_Ind_Vect_p);
    aug_nextE_Ind_Vect_p(nextlyInfE)=[];
    aug_Vect_Ls_p(nextlyInfE)=[];
    Lmin_offs_p(nextlyInfE)=[];

    %Store in cell array
    aug_nextE_Ind_Vect{p}=aug_nextE_Ind_Vect_p;
    aug_Vect_Ls{p}=aug_Vect_Ls_p;

    %Get index of subsequent next state that is off the grid
    x=1;

    %STEP 10: Create each P matrix
    %For P matrix, select rows corresponding to components in nextE_Ind_Vect
    %(Note: doing after P_fullmtx completed)
    for r=1:length(E_Ind_Vect_p)
        Ind_nextE=nextE_Ind_Vect_p(r);    %Get index of state stored in r-th row of nextE_Ind_Vect (i.e. the next energy state)

        %Get column number of next row of probabilities as RELATED to the NEXT ENERGY STATE INDEX (mapping to deterministic component!!!)
        c=find(abs(aug_nextE_Ind_Vect_p-Ind_nextE)<epsilon2,1); %Get from position of FIRST Ind_nextE in AUG_nextE_Ind_Vect!!!!! (b/c same width as AUGMENTED VECTOR)

        %IF off grid... STEP F
        if round(Ind_nextE)~=Ind_nextE
            %Count number of non-zero load probabilities for next state Ind_nextE
            nnzProb_nextE=numLoads_OffGrd_p(x);
            %Transform probabilities to custom distribution
            prob_nextE=ProbDistr_v2(nnzProb_nextE,feasStatesArr_p(r,3),(RES_L*(aug_Vect_Ls_p(c:c+nnzProb_nextE-1)-MIN_LOAD)+1),1);

            x=x+1;
        else %Otherwise...
            %Count number of non-zero probabilities in associated E-state row of P_fullmtx (i.e. Ind_nextE)
            nnzProb_nextE=nnz(P_fullmtx(Ind_nextE,:));      %Should be equal to number of repeats in nextE_Ind_Vect
            %Transform probabilities to custom distribution
            prob_nextE=ProbDistr_v2(nnzProb_nextE,feasStatesArr_p(r,3),(RES_L*(aug_Vect_Ls_p(c:c+nnzProb_nextE-1)-MIN_LOAD)+1),1);
        end

        %Store subscript pairs and associated values in array P
        len=length(prob_nextE); 
        cols=c:(c+nnzProb_nextE-1);
        P=[P;r*ones(len,1),cols',prob_nextE'];
    end

    %Store in p-th PF matrix, as well as in own P_mtx
    PF{p}=sparse(P(:,1),P(:,2),P(:,3)); %Store as SPARSE MATRIX 
    P_mtx{p}=PF{p};
    %Reset matrices/vectors
    P=[];
    aug_nextE_Ind_Vect_p=[];
    aug_Vect_Ls_p=[];


    else %ELSE IF CONTROLS ARE INFEASIBLE for this value of "p"
        %Ignore p
        PF{p}=[];
        P_mtx{p}=PF{p};
    end
  end




  %STEP 11: Construct each F matrix
  for p=1:p_max
      if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...

      aug_nextE_Ind_Vect_p=aug_nextE_Ind_Vect{p};
      aug_Vect_Ls_p=aug_Vect_Ls{p};
      offGrdNxtE1E2_p=offGrdNxtE1E2{p};
      %Index COLUMN of F matrix by ROW number of E_Ind_VectALL
      row=1; %Row being checked in E_Ind_VectALL to start when starting on next E_Ind vector. Reset it

      %Go through next E-state index vector for current value of p...
      for r=1:length(aug_nextE_Ind_Vect_p)

          %IF next state is off grid...
          if round(aug_nextE_Ind_Vect_p(r))~=aug_nextE_Ind_Vect_p(r)

              %STEP G:

              %Get individual off-grid next state indices
              idx=find(abs(offGrdNxtE1E2_p(:,3)-aug_nextE_Ind_Vect_p(r))<epsilon2);
              idx=idx(1); %Take only first element if repeats in offGrdNxtE1E2_p ................. ASSUMING REPEATED OFF-GRID linear INDICES MAP TO SAME SET OF index pairs

              nextE1_Ind=offGrdNxtE1E2_p(idx,1);
              nextE2_Ind=offGrdNxtE1E2_p(idx,2);
              nextL=aug_Vect_Ls_p(r);

              for col=1:length(E_Ind_VectALL)
                  %Get individual current state indices
                  E2_Ind=remainder(E_Ind_VectALL(col),N2);
                  E1_Ind=(E_Ind_VectALL(col)-E2_Ind)/N2+1;
                  L=E_VectALL_Ls(col);

                  %INTERPOLATION
                  %Check if (nextE1,nextE2) is on edge of square
                  %If so, apply different interpolation
                  if abs(nextL-L)<epsilon
                      if nextE1_Ind==E1_Ind
                          if floor(nextE2_Ind)==E2_Ind
                              q=1-(nextE2_Ind-E2_Ind);
                          elseif ceil(nextE2_Ind)==E2_Ind
                              q=1-(E2_Ind-nextE2_Ind);
                          else 
                              q=0;
                          end
                      elseif nextE2_Ind==E2_Ind
                          if floor(nextE1_Ind)==E1_Ind
                              q=1-(nextE1_Ind-E1_Ind);
                          elseif ceil(nextE1_Ind)==E1_Ind
                              q=1-(E1_Ind-nextE1_Ind);
                          else 
                              q=0;
                          end
                      %If on neither edge...
                      else 
                          %Check to find 4  points closest to (nextE1,nextE2) off grid.... FIND (E1_Ind, E2_Ind)
                          %CASE 1: round E1 down, round E2 down
                          if floor(nextE1_Ind)==E1_Ind && floor(nextE2_Ind)==E2_Ind
                              q=(1-(nextE1_Ind-E1_Ind))*(1-(nextE2_Ind-E2_Ind));
                          %CASE 2: round E1 up, round E2 down
                          elseif ceil(nextE1_Ind)==E1_Ind && floor(nextE2_Ind)==E2_Ind
                              q=(1-(E1_Ind-nextE1_Ind))*(1-(nextE2_Ind-E2_Ind));
                          %CASE 3: round E1 down, round E2 up
                          elseif floor(nextE1_Ind)==E1_Ind && ceil(nextE2_Ind)==E2_Ind
                              q=(1-(nextE1_Ind-E1_Ind))*(1-(E2_Ind-nextE2_Ind));
                          %CASE 4: round E1 up, round E2 up
                          elseif ceil(nextE1_Ind)==E1_Ind && ceil(nextE2_Ind)==E2_Ind
                              q=(1-(E1_Ind-nextE1_Ind))*(1-(E2_Ind-nextE2_Ind));
                          else
                             q=0; %If this state on grid not used for interpolation (not corner of encompassing square)
                          end
                      end
                  else
                      q=0;
                  end


                  %Store subscript pairs and associated weightings in F_p
                  if q~=0
                    F_p=[F_p;r,col,q]; %Use states on grid for interpolation, with WEIGHTING q
                  end
              end

          else %Otherwise, if ON-GRID...

              %If next state is currently infeasible...
              if aug_nextE_Ind_Vect_p(r)<E_Ind_VectALL(row) %(i.e. NOT continuously increasing in augmented vector)
                 row=1; %Restart from beginning of E_Ind_VectALL to find the state <----- ASSUMING ONLY 1 distinct new currently infeasible state!
              end

              while(E_Ind_VectALL(row)~=aug_nextE_Ind_Vect_p(r)) %While not reached mapping column in F (ONLY 1 per row)...
                  row=row+1;    %Continue
              end

              %Store subscript pairs and associated 1's (feasible next) in array F_p 
              F_p=[F_p;r,row,1]; %Mark next state as feasible
              row=min(row+1,length(E_Ind_VectALL)); %Start at next column in F next time, saturating at maximum
              %^-------- Assuming continuously increasing in augmented vector (fixed above)
          end
      end

      if isempty(F_p)   %If empty, ignore
         F{p}=0;
      else      %IN MOST CASES... 
        F{p}=sparse(F_p(:,1),F_p(:,2),F_p(:,3),max(F_p(:,1)),length(E_Ind_VectALL)); %STORE AS SPARSE MATRIX
        %Also add extra zeros at end to ensure dimensions of F{p} and E_Ind_VectALL match
      end

      %Correct for infeasible neighbouring states for next state interpolation in F matrix
      for i=1:size(F{p},1)
          if sum(F{p}(i,:))~=1
              %Distribute missing weight evenly to remaining neighbouring states
              inds=find(F{p}(i,:));%Indices of remaining states
              F{p}(i,inds)=F{p}(i,inds)+(1-sum(F{p}(i,:)))/nnz(F{p}(i,:));
          end
      end

      F_p=[]; %Reset

      if isempty(PF{p}) %If no next state..
          PF{p}=0;  %Ignore constraint
      else      %IN MOST CASES...
         PF{p} = PF{p}*F{p}; %Finish PF matrices 
      end



      else %ELSE IF CONTROLS ARE INFEASIBLE
        %Ignore p

        PF{p} = zeros(length(E_Ind_VectALL),length(E_Ind_VectALL));
      end
  end

  %STEP 12: Construct each G matrix
  for p=1:p_max
      if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...

      E_Ind_Vect_p=E_Ind_Vect{p};
      %Lmin_offs_p=Lmin_offs{p};
      E_Ind_L_Vect_p=E_Ind_L_Vect{p};
      feasStatesArr_p=feasStatesArr_ctrl{p};
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
          while(E_Ind_VectALL(round(row))~=E_Ind_Vect_p(r)) %While not reached mapping column in G (ONLY 1 per row)...
              row=row+1;    %Continue
          end
          %E_Ind_VectALL(row) is value (AND INDEX) for UNIQUE E-states
          if(boolNewEState==1)  %Only if distinct new state...
              row=row+1*( find(abs(E_VectALL_Ls(round(row):end)-E_Ind_L_Vect_p(r,2))<epsilon,1) - 1 ); %Add minimum load offset + negative loads offset to first state #
          else
              %Otherwise, do nothing because already starting from offset
          end

          %Store subscript pairs and associated 1's (feasible next) in array G_p 
          G_p=[G_p;r,row,1]; %Mark current state as feasible
          row=min(row+1,length(E_Ind_VectALL));  %Start at next column in G next time (continuously increasing)
      end

      if isempty(G_p)   %If empty, ignore constraint
         G{p}=0;
      else      %IN MOST CASES... 
         G{p}=sparse(G_p(:,1),round(G_p(:,2)),G_p(:,3),max(G_p(:,1)),length(E_Ind_VectALL)); %STORE AS SPARSE MATRIX
        %Also add extra zeros at end to ensure dimensions of G{p} and E_Ind_VectALL match
      end

      G_p=[]; %Reset


      else %ELSE IF CONTROLS ARE INFEASIBLE
        %Ignore p
        G{p} = zeros(length(E_Ind_VectALL),length(E_Ind_VectALL));
      end
  end

  %If get empty g vector...
  for p=1:p_max
    if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...

        if isempty(g{p})
            disp('ERROR!!'); %Error!!!
           g{p}=zeros(length(E_Ind_VectALL),1); %Set equal to zeros
        end

    else %If controls infeasible, set cost to infinity
        g{p}=Inf*ones(length(E_Ind_VectALL),1); %Set equal to infinity
    end
  end


  %% CORRECTED MATRICES (REMAINING infeasible states removed)
  %1) Coefficients

  Q=[];
  for p=1:p_max
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
  b=[]; b_dcpcvx=[];
  for p=1:p_max
      b=[b;g{p}];
      if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...
         b_dcpcvx=[b_dcpcvx;g{p}];
      else
         b_dcpcvx=[b_dcpcvx;zeros(length(E_Ind_VectALL),1)];
      end
  end

  %% PART B: COST APPROXIMATION AND BASIS FUNCTION GENERATION
  %Create initial design matrix (1 row per feasible state)
  %Use state aggregation for first 'num_steps', then start to generate
  %monomials

  feasE2s=[]; feasE1s=[]; feasLs=[]; %To add polynomial basis functions (order R-1)

    %{
    %1) find bounds of linear fit
    exprMax=0; exprMin=0;
    for i=1:N1
      E1=i-1;
      for j=1:N2
          E2=j-1;
          for k=1:size(feasStates,3)
                L=(k-1)/RES_L+MIN_LOAD;
                if(feasStates(i,j,k)==1)
                    exprMax=max(fitStateExpr(E1,E2,L),exprMax);
                    exprMin=min(fitStateExpr(E1,E2,L),exprMin);
                end
          end
      end
    end
    %2) determine suitable divisor so integer number steps
    num_steps=MAX_STEPS;
    while rem(exprMax-exprMin+1,num_steps) %While remainder, decrease steps
      num_steps=num_steps-1;
    end
    step_size=(exprMax-exprMin+1)/num_steps; %Get step (aggregation) size
    %3) finally, create feature vectors
    for i=1:N1
      E1=i-1;
      for j=1:N2
          E2=j-1;
          for k=1:size(feasStates,3)
                L=(k-1)/RES_L+MIN_LOAD;
                if(feasStates(i,j,k)==1)
                    %Create parameter fitting vector by aggregating states
                    %with same LINEAR EXPRESSION value (i.e. LINEAR FIT)
                    %Set expression
                    linStateExpr=fitStateExpr(E1,E2,L);
                    %Create vector
                    for l=1:num_steps
                        if ((exprMin+(l-1)*step_size)<=linStateExpr && ...
                                linStateExpr<(exprMin+l*step_size))
                            phi_vec=[zeros(1,l-1),1,zeros(1,num_steps-l)];
                        end
                    end
                    if phi_vec==0
                        disp('State out of bounds of aggregation!!');
                    end

                    %Add to design matrix, and then reset
                    Phi=[Phi;phi_vec];
                    phi_vec=0;


                    %PREPARE to form monomials for remainder of vectors
                    %If feasible state, add to list (for developing the
                    %design matrix of monomials)
                    feasE2s=[feasE2s;E2]; feasE1s=[feasE1s;E1]; feasLs=[feasLs;L];
                end
          end
      end
    end

    %Adjoin feasible state vectors to form a ?x3 array
    feasStatesArr=[feasE1s,feasE2s,feasLs];
    %
    %Create design matrix with fitting functions up to order R
    phi_poly=DesignMtx(feasStatesArr,ones(length(feasStatesArr),1),R);
    %Add to present basis vectors
    %OLD - use all vectors: Phi=phi_poly(:,1:end);
    lenPoly=size(phi_poly,2);
    %Determine the affine bases in order R polynomial by index
    indBaseE2=lenPoly-(R+2)+1; indBaseE1=lenPoly-(0.5*R^2+1.5*R+2)+1;
    %Get these columns of constant and linear terms
    boolNonAff=~[zeros(1,indBaseE1-1) 1 zeros(1,indBaseE2-indBaseE1-1) 1 zeros(1,lenPoly-2-indBaseE2) 1 1]; 
    %Remove these bases
    nonAff_phi_poly=phi_poly(:,boolNonAff);

    if rank(Phi)>1
        Phi=[Phi,nonAff_phi_poly]; %Ignore constant and linear terms (ALREADY ACCOUNTED FOR IN S.A.)
    else
        Phi=phi_poly;
    end

    origSizePhi=size(Phi,2); %Store old size of Phi
    %Remove linearly dependent column vectors UNTIL FULL RANK Phi
    oldPhi=Phi;
    Phi=licols(oldPhi,1e-13);

    %}
    %Alternative: use EXACT LP (Phi=I)
    %Phi=eye(length(E_Ind_VectALL));

    % Find state-relevance vector for minimization, c
    % TAKE c TO BE STEADY STATE ENTERING PROBABILITIES FOR EACH STATE
    % Probabilities are given in P_fullmtx (non-zero for feasible states)
    %trP_fullmtx=P_fullmtx';
    %c_state=trP_fullmtx(:); %Get probabilities for all states
    %c_state(c_state==0)=[]; %Remove zero probability states

    % Alternative: TAKE c TO BE ALL EQUAL WEIGHTS
    c_state=ones(length(E_Ind_VectALL),1);

    %Created LP matrices and vectors, so SOLVE OPTIMIZATION

    %cvx_solver Gurobi
    %Get approximate solution
    cvx_begin
    %cvx_solver_settings('Method',1) % Use dual simplex method
    cvx_solver_settings('Presolve',0) % Don't use presolver
    %cvx_solver_settings('FeasibilityTol',1e-4) %Set tolerance
    variable r_fit(length(E_Ind_VectALL))
    dual variables d
    %dual variables d2
    maximize( c_state'*r_fit ) %*Phi*r_fit ) % - gamma*norm(r_fit,1) )
    subject to
        %d : Q*Phi*r_fit <= b_dcpcvx
        d : Q*r_fit <= b_dcpcvx
        %d2 : Phi*r_fit >= 0
    cvx_end


    %    %Get error
    %    err=cost-Phi*r_fit;
    %    %Store NORMALIZED approximation error bases (2-NORM)
    %    approx_err=norm(err,2)/norm(cost,2);
    %    %approx_err=[approx_err;norm(cost-Phi*r_fit,2)/norm(cost,2)];
    %    %Store approximation
    %    approx=Phi*r_fit;
    %   
    %   %Plot approximate and actual costs
    %   figure
    %   plot(Phi*r_fit, '*'); hold on; plot(cost, '*');
    %   xlabel('State Index'); ylabel('Cost');
    %   legend('Approximate Cost','Actual Cost');
    %   title(strcat('Evaluating Approximation with',{' '},num2str(origSizePhi),'-Bases Fit, Rank of Phi=',num2str(rank(Phi))));

    %   %Get error in each state marginalized over E2 and L
    %   for r=1:length(feasStatesArr)
    %      E1=feasStatesArr(r,1)-E_MIN(1)+1;
    %      margApproxErr(E1)=margApproxErr(E1)+err(r);
    %   end
    %   %Plot
    %   figure
    %   plot(0:1:E_MAX(1),margApproxErr,'o');
    %   xlabel('Energy E1'); ylabel('Marginalized Error');
    %   title('Marginalized Error for default test');


    optD = d; %Get vector of FINAL dual
    %cost=Phi*r_fit; %Get FINAL approximated cost
    cost=r_fit; %Get FINAL exact cost


    %Format FINAL cost vector into E1xE2 matrices (one for each value of load)
    ConvCosts=FormatCostVect_v2(cost);

    %
    %% PART C: STATIONARY POLICY
    %PART 0: get state-action Q-values
    PF_mtx=[]; %Create aggregate transition matrix
    for p=1:p_max
        PF_mtx=[PF_mtx;PF{p}];
    end
    q_value=b+DISCOUNT*PF_mtx*cost; %DETERMINE OPTIMAL Q-VALUES


    %PART 1: Create vector of Q-values for each action
    %1) Augment q_value vector to include those of infeasible states too (infeasible q-value)
    %Make each E_Ind_Mtx same size to compare ALL states between different
    %values of p, once in vector form
    %->Augmented vectors for given values of p (i.e. like augmented versions of
    %E_Ind_VectALL, and subsets of E_MtxALL_Vect)
    E_MtxALL_Vect_subs={};
    maxNegL=size(E_Ind_Mtx{1},2); %Get maximum number of negative loads. <-------------------------------------------------------------------------------------------------------------------ASSUMING THIS CORRESPONDS TO p=1, since most charging controls
    for p=1:p_max
            E_Ind_Mtx_p=E_Ind_Mtx{p};
            E_Ind_Mtx_p=[zeros(size(E_Ind_Mtx_p,1),maxNegL-size(E_Ind_Mtx_p,2)),E_Ind_Mtx_p]; %Pad with zeros to indicate infeasibility of loads even BELOW MINLoad_E_state
            E_Ind_Mtx_p(:,size(E_Ind_Mtx_p,2)+1:size(E_Ind_MtxALL,2))=0; %Pad with zeros on side to make same size
            E_Ind_Mtx{p}=E_Ind_Mtx_p; %Overwrite with padded matrix
            %Convert to vector
            trE_Ind_Mtx_p=E_Ind_Mtx_p';
            E_MtxALL_Vect_subs{p}=trE_Ind_Mtx_p(:);
    end

    %2) Create augmented vectors of q-values for ALL states - feasible
    %AND INFEASIBLE TOO - for EACH CONTROL p

    %For each element in E_MtxALL_Vect_subs{p}, if...
    %a) 0, append INF to aug_optQ_subP{p} (since infeasible, for Q-value minimization)
    %b) non-zero, append some value from q_value to aug_optQ_subP{p}
    %where some value is next value in for i=1:p-1 sumLen=sumLen+len(vecti) end q_value(sumLen:sumLen+len(vectp))
    aug_optQ_subP={};
    indOptQ=1;
    for p=1:p_max
        aug_optQ_subP_p=[];
            E_MtxALL_Vect_subs_p=E_MtxALL_Vect_subs{p};
           for i=1:length(E_MtxALL_Vect_subs_p)
               if(E_MtxALL_Vect_subs_p(i)==0)
                  aug_optQ_subP_p=[aug_optQ_subP_p;INF_Q];
               else
                   %Find value in subvector of q-value just by continuously
                   %indexing through q-value in order <--------------------------Assuming optQ linearly indexed in order (E2, E1, L, D2, D1)
                   aug_optQ_subP_p=[aug_optQ_subP_p;q_value(indOptQ)];
                   indOptQ=indOptQ+1;
               end
           end
           aug_optQ_subP{p}=aug_optQ_subP_p;
    end

    %PART 2: Get optimal q-values vector
    %Create augmented q-value vector, for ALL states
    aug_optQ=[];
    for p=1:p_max
        if (~isempty(fbleEStates{p})) %IF CONTROLS NOT IMMEDIATELY INFEASIBLE for all states...
            aug_optQ=[aug_optQ;aug_optQ_subP{p}];
        else
            aug_optQ=[aug_optQ;INF_Q*ones(length(aug_optQ_subP{1}),1)];
        end
    end

    %Get final q-values for ALL states (augmented)
    aug_Q=aug_optQ;

    %Create augmented vector of all E_MtxALL_Vect_subs vectors
    aug_E_MtxALL_Vect=[];
    for p=1:p_max
        aug_E_MtxALL_Vect=[aug_E_MtxALL_Vect;E_MtxALL_Vect_subs{p}];
    end

    %Get FINAL q-values ONLY for feasible states (non-zero in aug_E_MtxALL_Vect)
    optQ=aug_Q(aug_E_MtxALL_Vect~=0);
    aug_optQ=aug_optQ(aug_E_MtxALL_Vect~=0);
  %}