boolSingleSequence=1;

if(~boolSingleSequence)
    %Visualize all possible policies
    E_Ind1=1:(E_MAX(1)-E_MIN(1)+1);
    E_Ind2=1:(E_MAX(2)-E_MIN(2)+1);
    indL=1:(MAX_LOAD-MIN_LOAD+1);

    ti=0;tf=LAST_ITER-1;

    figure
    %%% Layered-Surfaces Visualization %%%%
    for t=ti:tf
        subplot(1,LAST_ITER,t+1)
        for ind=1:length(indL)
            cost=V(:,:,ind,t+1);
            Z=(ind-1)*ones(length(E_Ind1),length(E_Ind2));
            surf(E_Ind2-1,E_Ind1-1,Z,cost)
            hold on
        end
        colorbar
        xlabel('E2');ylabel('E1');zlabel('L');
        title(['Value function at t=',num2str(t)]);
    end

    hold off;
    figure
    %%% 3D Scatter Plot w/ Labels %%%
    d=0.2; %displacement so the text does not overlay the data points
    for t=ti:tf
        subplot(1,LAST_ITER,t+1)
        for ind=1:length(indL)
            [X,Y,Z]=Cuboid(E_Ind2-1, E_Ind1-1,ind-1);
            scatter3(X,Y,Z)
            numPts=length(E_Ind2)*length(E_Ind1)*1;
            cost=reshape(V(:,:,ind,t+1),[numPts 1]);
            text(X,Y,Z+d,cellstr(num2str(cost)));
            hold on
        end
        xlabel('E2');ylabel('E1');zlabel('L');
        title(['Value function at t=',num2str(t)]);
    end

    hold off;


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
    
else
    %Visualize a possible sequence of loads, as previously found
    figure
    hold on;
    plot(optE1,'*','MarkerSize',10); plot(optE2,'*','MarkerSize',10); plot(Load,'o','MarkerSize',10);
    hold off;
    xlabel('Time');
    ylabel('Unit of energy');
    title('Change in energy stored vs. load');
    legend('E1','E2','Load');
    
    figure
    hold on;
    plot(D1Opt,'*','MarkerSize',10); plot(D2Opt,'*','MarkerSize',10); plot(Load,'o','MarkerSize',10);
    hold off;
    xlabel('Time');
    ylabel('Unit of energy');
    title('Optimal policy vs. load');
    legend('D1','D2','Load');
end