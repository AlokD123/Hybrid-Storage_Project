function [P_mtx,SOC_seq_disc,L_seq_disc,ctrs] = getRealProbMtx(SOC_seq,L_seq,numBinsE,numBinsL) %numBinsE should be N1*N2
%getRealProbMtx Uses real driving data, discretized into bins, to create demand probability distribution.
%Inputs:
%   SOC_seq: real battery SOC data. Assuming horizontal vector
%   L_seq: real energy demand data. Assuming horizontal vector
%   numBinsE: discretization of E1,E2 aggregate
%   numBinsL: discretization of L
%Outputs
%   P_mtx: matrix of probabilities for each state. Rows: E-states. Cols:  L-values

%1) Adjoin transposed sequences as column vectors (col1: SOC, col2: L)
%X=[SOC_seq',L_seq'];

%2) Get 2D histogram data
%[hist2D,ctrs]=hist3(X,'Nbins',[numBinsE,numBinsL]);

%
[~,SOC_edges]=discretize(SOC_seq,numBinsE);
[~,L_edges]=discretize(L_seq,numBinsL);


SOC_values = SOC_edges(2:end);
L_values=0.5*L_edges(1:end-1)+0.5*L_edges(2:end);
SOC_seq_disc=discretize(SOC_seq,SOC_edges,SOC_values);
L_seq_disc=discretize(L_seq,L_edges,L_values);
%}

%
L_seq_disc=L_seq_disc+max((SOC_seq_disc-L_seq_disc)-SOC_seq_disc(1),0);     %<------ ASSUMING SOC_seq_disc(1) is maximum value of E
L_seq_disc=L_seq_disc-max(SOC_seq_disc(end)-(SOC_seq_disc-L_seq_disc),0);   %<------ ASSUMING SOC_seq_disc(end) is minimum value of E

X=[SOC_seq_disc',L_seq_disc'];
[hist2D,ctrs]=hist3(X,'Nbins',[numBinsE,numBinsL]);

%{
[~,L_edges]=discretize(L_seq_disc,numBinsL);
L_values=L_edges(2:end);
L_seq_disc=discretize(L_seq_disc,L_edges,L_values);

%}


%Plot histogram
%hist3(X,'Nbins',[numBinsE,numBinsL]);

%3) Normalize over demands for each E-state
P_mtx=normalize(hist2D','norm',1)';

end