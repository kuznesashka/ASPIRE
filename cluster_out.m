function cluster_out  = cluster_out_spycirc(cluster, G3)

% -------------------------------------------------------------------------
% Reshape the cluster structure
% -------------------------------------------------------------------------
%   cluster -- structure [length(ind_m)x5], 
% first -- source location, 
% second -- spike timestamp, 
% third -- the subcorr value, 
% fourth -- index of the spike from the whole set (spike_ind)
% fifth -- spike_clust from Spyking Circus (number of template)
%
% OUTPUTS:
%   cluster -- structure 8 columns 
% first -- spike timestamp, 
% second -- Euclidean cluster, 
% third -- subcorr values from RAP-MUSIC, 
% fourth-sixth -- three coordinates of the spatial location of the dipole
% seventh -- vertex number
% eighth -- spike_clust from Spyking Circus (number of the template)
% 

cluster_out = [];
L = length(cluster);
for pre = 1:L

    stamps_pre = cluster{1,pre}(2,:)';
    gof        = cluster{1,pre}(3,:)';
    pos        = G3.GridLoc(cluster{1,pre}(1,:),:);
    vertex     = cluster{1,pre}(1,:)';
    SPC_cluster     = cluster{1,pre}(5,:)';

    cluster_out = [cluster_out; stamps_pre ones(size(stamps_pre))*pre gof pos vertex SPC_cluster];
end
