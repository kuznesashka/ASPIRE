function main_one_subject_propagation_run_from_python(cortex, G3, MRI, channels, paths, parameters)

% -------------------------------------------------------------------------
% All steps, one case
% -------------------------------------------------------------------------
% INPUTS:
%
% PATHS
%
%
% PARAMETERS
%
%
% OUTPUTS:
%
% sources - .mat file with sources
% clusters - .csv file with clusters
% results - .mat
% roc - .xlsx file with ROC curves
%
% _______________________________________________________
%

switch parameters.channel_type % channels you want to analyse ('grad' or 'mag')
    case 1, channel_type = 'grad';
        channel_idx     = setdiff(1:306, 3:3:306);
    case 2, channel_type = 'mag';
        channel_idx     = 3:3:306;
        
end
block_size = 600000; %ms
t1_IndMax = [];
t1_ValMax = [];
t1_ind_m = [];
t1_spikeind = [];
t1_spike_clust = [];
t1_spike_ind = [];

t3_IndMax = [];
t3_ValMax = [];
t3_ind_m = [];
t3_spikeind = [];
t3_spike_clust = [];
t3_spike_ind = [];

for data_n = 1:parameters.N_data
    switch data_n
        case 1, Data = load(parameters.Data_0);
        case 2, Data = load(parameters.Data_1);
        case 3, Data = load(parameters.Data_2);
    end
    block_begin = 41+block_size*(data_n-1);
    block_end   = block_size*(data_n)-30;


    [t1_IndMax, t1_ValMax, t1_spikeind, t1_spike_ind, t1_ind_m, t1_spike_clust] = ...
        dip_fit_one_block(t1_IndMax, t1_ValMax, t1_spikeind, t1_spike_ind, t1_ind_m, t1_spike_clust, ...
            block_size, block_begin, block_end, Data, G3, channels, channel_idx, t1, t2, t3, MRI, cortex)
    [t3_IndMax, t3_ValMax, t3_spikeind, t3_spike_ind, t3_ind_m, t3_spike_clust] = ...
        dip_fit_one_block(t3_IndMax, t3_ValMax, t3_spikeind, t3_spike_ind, t3_ind_m, t3_spike_clust, ...
            block_size, block_begin, block_end, Data, G3, channels, channel_idx, t1, t2, t3, MRI, cortex)

end

%% Clustering t1
IndMax      = t1_IndMax  
ValMax      = t1_ValMax  
ind_m       = t1_ind_m   
spikeind    = t1_spikeind
spike_clust = t1_spike_clust
spike_ind   = t1_spike_ind
save(paths.sources_saving_path_t1, 'IndMax','ValMax','ind_m','spikeind')
try
    cluster_t1 = clustering(spike_ind, ... 
                         G3, ...
                         parameters.clustering.N_MIN, ... 
                         ValMax, ...
                         IndMax, ...
                         ind_m, ...
                         parameters.clustering.THR_DIST, ... 
                         0, ... 
                         cortex, ... 
                         parameters.rap_music.RAP, ... 
                         spikeind, ... 
                         spike_clust);
catch
    if size(IndMax) ~= size(spike_ind)
        spike_ind  = spike_ind';
    end
    cluster_t1{1,1} = [IndMax(ind_m); spike_ind(ind_m); ValMax(ind_m)];
end

%% Clustering t3
IndMax      = t3_IndMax  
ValMax      = t3_ValMax  
ind_m       = t3_ind_m   
spikeind    = t3_spikeind
spike_clust = t3_spike_clust
spike_ind   = t3_spike_ind
save(paths.sources_saving_path_t3, 'IndMax','ValMax','ind_m','spikeind')
try
    cluster_t3 = clustering(spike_ind, ... 
                         G3, ...
                         parameters.clustering.N_MIN, ... 
                         ValMax, ...
                         IndMax, ...
                         ind_m, ...
                         parameters.clustering.THR_DIST, ... 
                         0, ... 
                         cortex, ... 
                         parameters.rap_music.RAP, ... 
                         spikeind, ... 
                         spike_clust);
catch
    if size(IndMax) ~= size(spike_ind)
        spike_ind  = spike_ind';
    end
    cluster_t3{1,1} = [IndMax(ind_m); spike_ind(ind_m); ValMax(ind_m)];
end

save(paths.clusters_saving_path ,'cluster')
save(paths.parameters_saving_path ,'parameters')
end

function [IndMax, ValMax, spikeind, spike_ind, ind_m, spike_clust] = ...
    dip_fit_one_block(IndMax, ValMax, spikeind, spike_ind, ind_m, spike_clust, ...
     block_size, block_begin, block_end, Data, G3, channels, channel_idx, t1, t2, t3, MRI, cortex)

    spc_data            = load(paths.detections);
    spike_ind_n         = int64(spc_data.spikes.ind');
    spike_clust_n       = spc_data.spikes.clusters';
    spike_clust_n       = spike_clust_n(spike_ind_n>block_begin & spike_ind_n<block_end);
    spike_ind_n         = spike_ind_n(spike_ind_n>block_begin & spike_ind_n<block_end) - block_size*(data_n-1);

    if isempty(spike_ind_n) ~= 1    
        [IndMax_n, ValMax_n, ind_m_n, spikeind_n] = dip_fit(spike_ind_n, Data, ...
            G3, channels, channel_idx, t1, t2, t3, MRI, cortex);

        l = length(spikeind);
        IndMax   = [IndMax IndMax_n];
        ValMax   = [ValMax ValMax_n];
        spikeind = [spikeind; spikeind_n    + block_size*(data_n-1)];
        spike_ind = [spike_ind; spike_ind_n + block_size*(data_n-1)];
        ind_m    = [ind_m  ind_m_n  + l(1)];
        spike_clust = [spike_clust; spike_clust_n];
    end
end

