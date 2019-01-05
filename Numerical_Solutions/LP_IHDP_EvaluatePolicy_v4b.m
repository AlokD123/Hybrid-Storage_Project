%GET INFINITE HORIZON POLICY USING LP SOLUTION
%V4: COMBINED CONTROLS

%NumIter=10; %Number of iterations of the policy to do

global MAX_CHARGE; global MAX_DISCHARGE;
global MAX_LOAD; global MIN_LOAD;

%E_MIN=[0;0]; %Minimum energy to be stored (lower bound)
%E_MAX=[5;4]; %Maximum energy to be stored (upper bound)
%E1_INIT=E_MAX(1); 
%E2_INIT=E_MAX(2);

%%Evaluate policy for a random sequence of loads (ONLINE)
%Set up matrices
optE1(1)=E1_INIT; optE2(1)=E2_INIT;
% Debug counts
countOOB=0;         %Out of bounds count
countRepeatZeros=0; %Count of repeated zero loads

%LOAD SEQUENCE!!!!
seqL=[]; %Reset load sequence

%FLUCTUATING LOAD
%
seqL=[0.25 1 0.75 0 0.5 1 1.75 2 1.25 2 2.5 2.75 1.5 3 2.25 2.75 3.5 3 3.25 3.25]; %FLUCTUATING LOAD
for i=1:1
    seqL=-[0.5*(1+normrnd(0.4,0.4,1,140)).*(1+0.2*sin(linspace(1,50,140))).*(-0.5*sin(linspace(1,10,140)))];

end
%}

%CYCLIC LOAD
%{
seqL=[0.25 1 0.75 0 0.5 1 1.75 2 1.25 2 2.5 2.75 1.5 3 2.25 2.75 3.5 3 3.25 3.25]; %FLUCTUATING LOAD
for i=1:1
    seqL=[seqL,fliplr(seqL),-seqL,-fliplr(seqL),seqL,fliplr(seqL),-seqL,-fliplr(seqL)];

end
%}

%SMOOTH LOAD
%{
seqL=[0.25 0 0.25 0.5 0.25 0.75 1.25 0.5 1.75 2 2.25 1.75 2.25 2.75 2.75 3 2.75 3 3.25 3.5];
for i=1:1
    seqL=[seqL,fliplr(seqL(end/2:end)).*linspace(1,1.75,length(seqL(end/2:end)))];
    seqL=[seqL,(1+normrnd(0.2,0.1,1,129)).*linspace(3.2,-3,129).*(1+0.5*sin(linspace(1,10,129)))];
end
%}

%RAMP LOAD
%{
seqL=[0.25 0 0.25 0.5 0.75 1 1 0.75 1 1.25 1.25 1.25 1 1.5 1.5 1.75 2 2.25 2.25 2.25];
for i=1:1
    seqL=[seqL,2*rand(1)*seqL(end/2:end)+0.5];
    seqL=[seqL,rand(1)*fliplr(seqL(end/2:end))];
    seqL=[seqL/5,-rand(1)*seqL(end/2:end)/2];
    seqL=[seqL,-rand(1)*fliplr(seqL(1:end))/2];
    seqL=[seqL,fliplr(seqL(end-17:end))+0.25]*1;
end
%}



%seqL=[testLoad,0];

seqL=seqL/1;
seqL(seqL>MAX_LOAD)=MAX_LOAD;
seqL(seqL<MIN_LOAD)=MIN_LOAD;


plot(seqL);
ylim([-5 8])
%}

NumIter=length(seqL); %Number of iterations of the policy to do

DeltaL_min=0.2;

t_ind_VI=1; %Start evaluation
while t_ind_VI<NumIter
    %Set state index
    indE1=optE1(t_ind_VI)-E_MIN(1)+1;
    indE2=optE2(t_ind_VI)-E_MIN(2)+1;
    
    %% Create demand sequence
    %OPTION 1: UNUSED
    %Create random demand from IID Uniform probability sequence
    MAX_LOAD_STATE=optE1(t_ind_VI)+optE2(t_ind_VI)-1; %Maximum possible load limited to total energy stored in that state
    if(MAX_LOAD_STATE==Inf)
        MAX_LOAD_STATE=MAX_LOAD;
    end
    %randL=randi(MAX_LOAD_STATE-MIN_LOAD+1,1,1)+MIN_LOAD-1;

    %Option 2: create sample sequence of pseudo-random demands
    %Select between sequences
    %L=randL;
    L= seqL(t_ind_VI);
    %L=min(t_ind_VI,MAX_LOAD_STATE); %Ramp
   
    %% Run sequences online
    
    %MOST IMPORTANT: RUN POLICY ONLINE!
    %Get optimal controls for given state
    U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
    
    %While load makes for infeasible state, change till feasible. L=0
    %always feasible, so break if get there.
    while isempty(U1)||(L-U1)<-MAX_CHARGE(2)-eps||(L-U1)>MAX_DISCHARGE(2)+eps %Infeasible if no optimal control OR C2 constraints violated
        if L<0
            L=L+DeltaL_min;
            %L=0;                   %To test with no R.B.
        elseif L>0
            L=L-DeltaL_min;
            %L=max(0,L-DeltaL_min); %To test with no R.B.
        else
            U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
            break;
        end
        U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
    end
    
    
    %Calculate the state these values of u and w will lead to, even if impossible...
    [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
        
    %If leads to next state guaranteed out of bounds, change load
    while nextE1>E_MAX(1)+eps||nextE1<E_MIN(1)-eps||nextE2>E_MAX(2)+eps||nextE2<E_MIN(2)-eps||(L-U1)<-MAX_CHARGE(2)-eps||(L-U1)>MAX_DISCHARGE(2)+eps
        %Increment out of bounds count
        countOOB=countOOB+1;
        %Change load till feasible
        if L<0
            L=L+DeltaL_min;
            %L=0;                   %To test with no R.B.
        elseif L>0
            L=L-DeltaL_min;
            %L=max(0,L-DeltaL_min); %To test with no R.B.
        else            %For L=0, just get controls and exit
            U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
            [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
            break;
        end
        
        %Calculate next state for new load value
        U1=GetPOpt_wRegenB_v3(indE1,indE2,L);
        [nextE1,nextE2]=optNextStateLimited_v3(optE1(t_ind_VI),optE2(t_ind_VI),U1,L);
    end
    if(L==0 && countOOB~=0)
        countRepeatZeros=countRepeatZeros+1;
    end
    Load(t_ind_VI)=L;          %Hold value of load (for reference)

    
    %Get specific control sequence for given load sequence
    U1Opt(t_ind_VI)=U1;
    %Get optimal next state, starting at E_INIT
    optE1(t_ind_VI+1)=nextE1;
    optE2(t_ind_VI+1)=nextE2;
    %Note: cross-check with optNextStateLtd?
    
    %If counted zero loads too many times, break loop
    if(countRepeatZeros==MAX_NUM_ZEROS)
        optE1(t_ind_VI+1)=optE1(t_ind_VI); optE2(t_ind_VI+1)=optE2(t_ind_VI);
        Load(t_ind_VI+1)=0;
    end
    
    %Finish iteration
    fbleCtrlOnGrd=[];
    fprintf("Completed iteration t=%d\n",t_ind_VI);
    
    t_ind_VI=t_ind_VI+1;
end