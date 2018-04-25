clear ERRORS_LP;
INFCOST=max(max(max(ConvCosts)));
ERRORS_LP=Costs(:,:,:)-ConvCosts(:,:,:);
%ERRORS_LP(abs(ERRORS_LP)>(INFCOST-0.1))=0;
%ERRORS_LP=ERRORS_LP./ConvCosts;

ERRORS_LP(ERRORS_LP==Inf)=+100;
ERRORS_LP(ERRORS_LP~=+100)=0;

ERRORS_LP(end+1,:,:)=0;
ERRORS_LP(:,end+1,:)=0;

if(1)
    %Visualize all possible policies
    E_Ind1=1:(E_MAX(1)-E_MIN(1)+2);
    E_Ind2=1:(E_MAX(2)-E_MIN(2)+2);
    indL=1:(MAX_LOAD-MIN_LOAD+1);

    %ti=0;tf=MAX_ITER-1;

    figure
    %%% Layered-Surfaces Visualization %%%%
    for t=1:1
        %subplot(1,1,t+1)
        for ind=1:length(indL)
            abserror=ERRORS_LP(:,:,ind);
            Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
            surf(E_Ind2-1,E_Ind1-1,Z,abserror)
            hold on
        end
        colorbar
        title(colorbar,'Error (+100 indicates infeasible)')
        xlabel('Supercapacitor Energy (E2)');ylabel('Battery Energy (E1)');zlabel('Demand (L)');
        title(['Percent Error between Costs in IHDP and Reduced LP']);
    end

%     hold off;
%     figure
%     %%% 3D Scatter Plot w/ Labels %%%
%     d=0.2; %displacement so the text does not overlay the data points
%     for t=ti:tf
%         subplot(1,LAST_ITER,t+1)
%         for ind=1:length(indL)
%             [X,Y,Z]=Cuboid(E_Ind2-1, E_Ind1-1,ind-1);
%             scatter3(X,Y,Z)
%             numPts=length(E_Ind2)*length(E_Ind1)*1;
%             cost=reshape(V(:,:,ind,t+1),[numPts 1]);
%             text(X,Y,Z+d,cellstr(num2str(cost)));
%             hold on
%         end
%         xlabel('E2');ylabel('E1');zlabel('L');
%         title(['Value function at t=',num2str(t)]);
%     end

%    hold off;


    %%%% NOT WORKING.... Array-of-Graphs Visualization %%%%
    % len1=length(E_Ind2);
    % len2=length(indL);
    %Plot E2 vs Cost
    % for x1=1:len1
    %     for x2=1:len2
    %         subplot(len1,len2, sub2ind([len1,len2],x1,x2));
    %         plot(E_Ind2,V(x1,:,x2,1))
    %         axis([1 length(E_Ind1)+1 0 10])
    %         title(['E1=',num2str(x1-1),', L=',num2str(x2-1)]);
    %     end
    % end


    %%% INCOMPLETE ... 3D Scatter Plot w/ Colour - Visualization %%%
    % X=E_Ind2-1;
    % for ind=1:length(E_Ind1)
    %     X=[X X];
    % end
    % %%%% INCOMPLETE
    % for ind=1:1
    %     %cost=V(:,:,ind,1);
    %     scatter3(E_Ind2-1,E_Ind1-1,)
    %     hold on
    % end
    % xlabel('E2');ylabel('E1');zlabel('L');
    % title('Cost');
    
end