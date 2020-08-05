%Script to visualize cost of applied controls at a point in time
%UNCOMBINED controls (non-negatives)

clearvars -except seqL MAX_DISCHARGE ALPHA_D ALPHA_C K

D1vals=1:MAX_DISCHARGE(1);
D2vals=1:MAX_DISCHARGE(2);
Lvals=1:(MAX_DISCHARGE(1)+MAX_DISCHARGE(2));

for i=1:(MAX_DISCHARGE(1))
    for j=1:(MAX_DISCHARGE(2))
        for l=1:(MAX_DISCHARGE(1)+MAX_DISCHARGE(2))
            ctrlCost1(i,j,l)=K*(i).^2+((1-ALPHA_D(1))*(i)+(1-ALPHA_D(2))*(j)+(1-ALPHA_C(2))*(i+j-l));
        end
    end
end

figure
for ind=1:(length(Lvals))
    cost=ctrlCost1(:,:,ind);
    Z=(ind-1)*ones(length(D1vals),length(D2vals));
    surf(D2vals,D1vals,Z,cost)
    hold on
end
colorbar
xlabel('D2');ylabel('D1');zlabel('L');
title(['Control cost function']);

figure
%%% 3D Scatter Plot w/ Labels %%%
d=0.2; %displacement so the text does not overlay the data points
for ind=1:length(Lvals)
    [X,Y,Z]=Cuboid(D2vals, D1vals,ind-1);
    scatter3(X,Y,Z)
    numPts=(length(D2vals))*(length(D1vals))*1;
    cost=reshape(ctrlCost1(:,:,ind),[numPts 1]);
    text(X,Y,Z+d,cellstr(num2str(cost)));
    hold on
end
xlabel('D2');ylabel('D1');zlabel('L');
title(['Control cost function']);