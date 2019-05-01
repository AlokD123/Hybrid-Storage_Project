function [P_fullmtx,SOC_seq_disc,L_seq_disc,ctrs] = getRealProbMtx(SOC_seq,L_seq,binsE,E_Ind_MtxALL) %numBinsE should be N1*N2
%getRealProbMtx Uses real driving data, discretized into bins, to create demand probability distribution.
%Inputs:
%   SOC_seq: real battery energy data. Assuming horizontal vector
%   L_seq: real energy demand data. Assuming horizontal vector
%   binsE: discretization of total energy (centre values of bins)
%   E_Ind_MtxALL: matrix indicating if state is feasible
%Outputs:
%   P_fullmtx: matrix of probabilities for each state. Rows: E-states. Cols:  L-values

global RES_L; global SCALE_BATT;

%1) Adjoin transposed sequences as column vectors (col1: SOC, col2: L)
%X=[SOC_seq',L_seq'];

%2) Get 2D histogram data
%[hist2D,ctrs]=hist3(X,'Nbins',[numBinsE,numBinsL]);

%
SOC_edges = conv(binsE, [0.5, 0.5], 'valid');
SOC_edges = [binsE(1)+(binsE(1)-SOC_edges(1));SOC_edges;binsE(end)+(binsE(end)-SOC_edges(end))];

L_max_abs=max(max(L_seq),-min(L_seq));
L_edges=(floor(min(-L_max_abs,L_max_abs)*RES_L)/RES_L-0.5/RES_L) : 1/RES_L: (ceil(max(-L_max_abs,L_max_abs)*RES_L)/RES_L+0.5/RES_L);

SOC_values = binsE;
L_values=0.5*L_edges(1:end-1)+0.5*L_edges(2:end);
SOC_seq_disc=discretize(SOC_seq,SOC_edges,SOC_values);
L_seq_disc=discretize(L_seq,L_edges,L_values);
%}

[histbefore,~]=hist3([SOC_seq_disc',L_seq_disc'],'Ctrs',{SOC_values L_values});
%
L_seq_disc=L_seq_disc+max((SOC_seq_disc-L_seq_disc)-SOC_seq_disc(1),0);     %<------ ASSUMING SOC_seq_disc(1) is maximum value of E
L_seq_disc=L_seq_disc-max(SOC_seq_disc(end)-(SOC_seq_disc-L_seq_disc),0);   %<------ ASSUMING SOC_seq_disc(end) is minimum value of E

X=[SOC_seq_disc',L_seq_disc'];
%[hist2D,ctrs]=hist3(X,'Nbins',[numBinsE,numBinsL]);
[hist2D,ctrs]=hist3(X,'Ctrs',{SOC_values L_values});

fprintf('Modified %f percent of data\n\n',(1-nnz(hist2D)/nnz(histbefore))*100);

%hist3(X,'Ctrs',{SOC_values L_values}); set(gca,'XTickLabel',linspace(SCALE_BATT*0.31,SCALE_BATT*0.76,9)); xlabel('Total Energy Stored (kWh)'); ylabel('Demand (Wh)'); zlabel('Count'); title('Histogram of stored energy and demand');



%Add zero prob states
missingLoadsNum=size(E_Ind_MtxALL,2)-size(hist2D,2);
hist2D=[zeros(size(E_Ind_MtxALL,1),missingLoadsNum/2),hist2D,zeros(size(E_Ind_MtxALL,1),missingLoadsNum/2)];

%For inappropriate demands in data (based on this discretization), mask out
fprintf('Lost %f percent of data\n\n',sum(sum(hist2D.*~E_Ind_MtxALL))/sum(sum(hist2D))*100);
hist2D=hist2D.*E_Ind_MtxALL;

%{
hist=load('hist4.mat');
hist2D=hist.hist2D;
%}

%3) Normalize over demands for each E-state
P_fullmtx=normalize(hist2D','norm',1)';



end