function [L_seq,Tot_E_seq] = useDrivingData(P_seq,E_norm_seq)
%useDrivingData Converts normalized energy to unitless energy and power to unitless demand, for use in driving problem.

global E_MAX; global MAX_LOAD; global MIN_LOAD;

L_seq=P_seq/max(max(P_seq),-min(P_seq))*max(MAX_LOAD,-MIN_LOAD);
Tot_E_seq=(E_norm_seq-min(E_norm_seq))/max(E_norm_seq-min(E_norm_seq))*sum(E_MAX);  %Total energy (unitless, scaled)

end