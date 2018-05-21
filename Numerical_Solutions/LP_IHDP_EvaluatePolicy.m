%GET INFINITE HORIZON POLICY USING LP SOLUTION
NumIter=10; %Number of iterations of the policy to do

%E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
%E_MAX=[5;4]; %Maximum energy to be stored (upper bound)
%E1_INIT=E_MAX(1); 
%E2_INIT=E_MAX(2);

%%STEP 2: evaluate policy for a random sequence of loads (ONLINE)
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads

%LOAD SEQUENCE!!!!
seqL=[]; %Reset load sequence
seqL=ones(NumIter,1); %CONSTANT LOAD
%seqL=[3 1 0 1 0 2 0 1 1 0]; %FLUCTUATING LOAD
%seqL=[3 0 1 0 1 0 2 0 1 0];

t_ind_VI=1; %Start evaluation
while t_ind_VI<NumIter
    %Set state index
    indE1=optE1(t_ind_VI)-E_MIN(1)+1;
    indE2=optE2(t_ind_VI)-E_MIN(2)+1;
    
    %UNUSED
    %Create random demand from IID Uniform probability sequence
    MAX_LOAD_STATE=optE1(t_ind_VI)+optE2(t_ind_VI)-1; %Maximum possible load limited to total energy stored in that state
    if(MAX_LOAD_STATE==Inf)
        MAX_LOAD_STATE=MAX_LOAD;
    end
    randL=randi(MAX_LOAD_STATE-MIN_LOAD+1,1,1)+MIN_LOAD-1;

    %Or, use sample sequence of pseudo-random demands
    %Select between sequences
    %L=randL;
    L= seqL(t_ind_VI);
    %L=min(t_ind_VI,MAX_LOAD_STATE); %Ramp
    
    %While load makes for infeasible state, decrease till feasible
    while(isempty(find(fullPolicyMtx(indE1,indE2,L-MIN_LOAD+1,:)==1,1)))
        L=L-1;
        if(L<0)
            L=0;
            break;
        end
    end
    %Get optimal p for given state
    pOpt=find(fullPolicyMtx(indE1,indE2,L-MIN_LOAD+1,:)==1);
    %Get associated optimal controls
    D2_Ind=remainder(pOpt,P2);
    D1_Ind=(pOpt-D2_Ind)/P2+1;
    D1=D1_Ind-1; D2=D2_Ind-1;
    
    %Calculate the state these values of u and w will lead to, even if impossible...
    [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
        
    %If leads to next state guaranteed out of bounds, decrease load
    while(nextE1>E_MAX(1)||nextE1<E_MIN(1)||nextE2>E_MAX(2)||nextE2<E_MIN(2)||(D1+D2-L)<0)
        L=L-1;%Decrement load
        %Increment out of bounds count
        countOOB=countOOB+1;
        if(L<0)
            L=0;
            break;
        end
        
        %Calculate next state for new load value
        pOpt=find(fullPolicyMtx(indE1,indE2,L-MIN_LOAD+1,:)==1);
        D2_Ind=remainder(pOpt,P2);
        D1_Ind=(pOpt-D2_Ind)/P2+1;
        D1=D1_Ind-1; D2=D2_Ind-1;
        [nextE1,nextE2]=optNextStateLimited(optE1(t_ind_VI),optE2(t_ind_VI),D1,D2,L);
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    indL=L-MIN_LOAD+1;
    Load(t_ind_VI)=L;          %Hold value of load (for reference)
    %Get costs of other possible states for this state and load??
    %(For comparison)

    %Round next state to nearest int
    nextE1=round(nextE1);
    nextE2=round(nextE2);
    %TO DO: APPLY STATE INTERPOLATION for better accuracy
    
    %Get specific control sequence for given load sequence
    D1Opt(t_ind_VI)=D1;
    D2Opt(t_ind_VI)=D2;
    %Get optimal next state, starting at E_INIT
    optE1(t_ind_VI+1)=nextE1;
    optE2(t_ind_VI+1)=nextE2;
    %Note: cross-check with optNextStateLtd?
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        optE1(t_ind_VI+1)=optE1(t_ind_VI); optE2(t_ind_VI+1)=optE2(t_ind_VI);
        Load(t_ind_VI+1)=0;
    end
    t_ind_VI=t_ind_VI+1;
end