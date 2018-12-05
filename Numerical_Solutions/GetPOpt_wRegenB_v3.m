function [U1] = GetPOpt_wRegenB_v3(indCurrE1,indCurrE2,CurrL)
%GetPOpt Returns the optimal control for the given state, interpolating
%when the state is off the grid
%   Input: state
%   Output: optimal (battery) control

global epsilon2; global E_Ind_VectALL; global E_VectALL_Ls; global N2;
global fullPolicyMtx; global MIN_LOAD;

if (round(indCurrE1)~=indCurrE1)||(round(indCurrE2)~=indCurrE2)||(round(CurrL)~=CurrL) %If state falls off the grid...
    q_pOpt_Mtx=[];
    
   %A) Go through all states and...
   for i=1:length(E_Ind_VectALL)
      %1) Get individual current state indices
      E2_Ind=remainder(E_Ind_VectALL(i),N2);
      E1_Ind=(E_Ind_VectALL(i)-E2_Ind)/N2+1;
      L=E_VectALL_Ls(i);

      %2) GET INTERPOLATION WEIGHTINGS
      %Check if (nextE1,nextE2) is on edge of square
      %If so, apply different interpolation
      if CurrL==L
          if indCurrE1==E1_Ind
              if floor(indCurrE2)==E2_Ind
                  q=1-(indCurrE2-E2_Ind);
              elseif ceil(indCurrE2)==E2_Ind
                  q=1-(E2_Ind-indCurrE2);
              else 
                  q=0;
              end
          elseif indCurrE2==E2_Ind
              if floor(indCurrE1)==E1_Ind
                  q=1-(indCurrE1-E1_Ind);
              elseif ceil(indCurrE1)==E1_Ind
                  q=1-(E1_Ind-indCurrE1);
              else 
                  q=0;
              end
          %If on neither edge...
          else 
              %Check to find 4  points closest to (nextE1,nextE2) off grid.... FIND (E1_Ind, E2_Ind)
              %CASE 1: round E1 down, round E2 down
              if floor(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind
                  q=(1-(indCurrE1-E1_Ind))*(1-(indCurrE2-E2_Ind));
              %CASE 2: round E1 up, round E2 down
              elseif ceil(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind
                  q=(1-(E1_Ind-indCurrE1))*(1-(indCurrE2-E2_Ind));
              %CASE 3: round E1 down, round E2 up
              elseif floor(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind
                  q=(1-(indCurrE1-E1_Ind))*(1-(E2_Ind-indCurrE2));
              %CASE 4: round E1 up, round E2 up
              elseif ceil(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind
                  q=(1-(E1_Ind-indCurrE1))*(1-(E2_Ind-indCurrE2));
              else
                 q=0; %If this state on grid not used for interpolation (not corner of encompassing square)
              end
          end
      else
        %Check to find 8 points closest to (CurrE1,CurrE2,CurrL) off grid.... FIND (E1_Ind, E2_Ind, L)
        %CASE 1: round E1 down, round E2 down, round L down
        if floor(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind && floor(CurrL)==L
          q=(1-(indCurrE1-E1_Ind))*(1-(indCurrE2-E2_Ind))*(1-(CurrL-L));
        %CASE 2: round E1 up, round E2 down, round L down
        elseif ceil(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind && floor(CurrL)==L
          q=(1-(E1_Ind-indCurrE1))*(1-(indCurrE2-E2_Ind))*(1-(CurrL-L));
        %CASE 3: round E1 down, round E2 up, round L down
        elseif floor(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind && floor(CurrL)==L
          q=(1-(indCurrE1-E1_Ind))*(1-(E2_Ind-indCurrE2))*(1-(CurrL-L));
        %CASE 4: round E1 up, round E2 up, round L down
        elseif ceil(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind && floor(CurrL)==L
          q=(1-(E1_Ind-indCurrE1))*(1-(E2_Ind-indCurrE2))*(1-(CurrL-L));
        %CASE 5: round E1 down, round E2 down, round L up
        elseif floor(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind && ceil(CurrL)==L
          q=(1-(indCurrE1-E1_Ind))*(1-(indCurrE2-E2_Ind))*(1-(L-CurrL));
        %CASE 6: round E1 up, round E2 down, round L up
        elseif ceil(indCurrE1)==E1_Ind && floor(indCurrE2)==E2_Ind && ceil(CurrL)==L
          q=(1-(E1_Ind-indCurrE1))*(1-(indCurrE2-E2_Ind))*(1-(L-CurrL));
        %CASE 7: round E1 down, round E2 up, round L up
        elseif floor(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind && ceil(CurrL)==L
          q=(1-(indCurrE1-E1_Ind))*(1-(E2_Ind-indCurrE2))*(1-(L-CurrL));
        %CASE 8: round E1 up, round E2 up, round L up
        elseif ceil(indCurrE1)==E1_Ind && ceil(indCurrE2)==E2_Ind && ceil(CurrL)==L
          q=(1-(E1_Ind-indCurrE1))*(1-(E2_Ind-indCurrE2))*(1-(L-CurrL));
        else
         q=0; %If this state on grid not used for interpolation (not corner of encompassing square)
        end
      end


      %3) Store states ON GRID and ASSOCIATED weightings in q_pOpt_Mtx
        q_pOpt_Mtx=[q_pOpt_Mtx;E1_Ind,E2_Ind,L,q]; %Use states on grid for interpolation, with WEIGHTING q
  end
    
   %B) Correct for infeasible neighbouring states for state interpolation in q_pOpt_Mtx
      if abs(sum(q_pOpt_Mtx(:,4))-1)>epsilon2
          %Get missing weight
          missingWeight=sum(q_pOpt_Mtx(:,4));
          %Distribute missing weight evenly to remaining neighbouring states
          for j=1:size(q_pOpt_Mtx,1)
              if q_pOpt_Mtx(j,4)~=0 %If neighbouring, i.e. non-zero associated weight ...
                 q_pOpt_Mtx(j,4)=q_pOpt_Mtx(j,4)+(1-missingWeight)/nnz(q_pOpt_Mtx(:,4)); %Add missing weight
              end
          end
      end
   
   %C) FINALLY, calculate optimal control value using interpolation
      U1_interp=0;
      for i=1:length(E_Ind_VectALL) %For every state on grid...
          if max(fullPolicyMtx(q_pOpt_Mtx(i,1),q_pOpt_Mtx(i,2),q_pOpt_Mtx(i,3)-MIN_LOAD+1,:))==0 %If demand is OUTSIDE the space of the grid (not just off of the grid)...
              U1_onGRD=[];
          else
            U1_onGRD=GetPOpt_wo_Interp_wRegenB_v3(q_pOpt_Mtx(i,1),q_pOpt_Mtx(i,2),q_pOpt_Mtx(i,3)); %Get optimal controls on grid
          end
          
        %If no optimal control value for state, due to approximation error..... <---------------------------------------------------------ERROR!!!!
        if isempty(U1_onGRD) %<--------------------------------------------------------------------------------------------------------Workaround!!
            if q_pOpt_Mtx(i,4)>0
                disp('Approx Error: U1!');
                break;
            end
            D1_onGRD=0;
        end
        %Create weighted sum of optimal controls at all points in grid
        U1_interp=U1_interp+q_pOpt_Mtx(i,4)*U1_onGRD;
      end
      
      U1=U1_interp(1); %Output interpolated optimal control
      
else %If state is on the grid
    u=GetPOpt_wo_Interp_wRegenB_v3(indCurrE1,indCurrE2,CurrL); %Get optimal control normally
    U1=u(1);
end



end

