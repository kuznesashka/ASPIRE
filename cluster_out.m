function cluster_out  = cluster_out_spycirc(cluster, G3)

% export in csv format 2 columns: time in ms, number of cluster
% cluster_out = 
% [stamps_pre cluster gof pos_xyz vertex]

cluster_out = []
L = length(cluster);
     for pre = 1:L
        stamps_pre = cluster{1,pre}(2,:)'; % column
        gof        = cluster{1,pre}(3,:)';
        pos        = G3.GridLoc(cluster{1,pre}(1,:),:);
        vertex     = cluster{1,pre}(1,:)';
        cluster_out = [cluster_out;
            stamps_pre ones(size(stamps_pre))*pre gof pos vertex];
     end
     