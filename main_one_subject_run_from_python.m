function main_one_subject_run_from_python(cortex,  Data, G3, MRI, channels, paths, parameters)

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
IndMax = [];
ValMax = [];
ind_m = [];
spikeind = [];
for data_n = 1:paths_params.N_data
    switch data_n
        case 1, Data = load(paths_params.Data_0);
        case 2, Data = load(paths_params.Data_1);
        case 3, Data = load(paths_params.Data_2);
    end
    block_begin = 41+block_size*(data_n-1);
    block_end   = block_size*(data_n)-30;
    %% 2. Spike detection
    switch parameters.spikes_detection
        case 1 % visual markings
            %manual_data = csvread(paths.path_vis_detections);
            manual_data       = load(paths.detections);
            picked_components = []; % ICA relevant
            picked_comp_top   = []; % ICA relevant
            spcirc_clust      = []; % SPC relevant (maybe delete)
            spike_ind         = manual_data(:,1);
            spike_ind         = spike_ind(spike_ind<600-30)*1000;

        case 2 % ICA based
            if parameters.run_ica
                [spike_ind, picked_components, picked_comp_top,component_indicatior] = ...
                    ICA_detection(Data, ...
                                  G3, ...
                                  channel_type, ...
                                  parameters.detection.ICA.decision, ...
                                  parameters.detection.ICA.f_low, ...
                                  parameters.detection.ICA.f_high);
                save(paths.detections,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
            end
            load(paths.detections,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
            spcirc_clust       = [];
            spike_clust        = zeros(size(spike_ind))';
            spike_clust        = spike_clust(spike_ind<600000-30 & spike_ind>41);
            spike_ind          = spike_ind(spike_ind<600000-30 & spike_ind>41);

        case 3 % Spiking circus based
            %spcirc_data = csvread([paths.path_SPC_detections,'_', channel_type, '.csv'],1,0);
            spc_data          = load(paths.detections);
            picked_components = []; % ICA relevant
            picked_comp_top   = []; % ICA relevant
            spike_ind         = spc_data.spikes.ind';
            spike_clust       = spc_data.spikes.clusters';

            spike_clust       = spike_clust(spike_ind>block_begin & spike_ind<block_end);
            spike_ind         = spike_ind(spike_ind>block_begin & spike_ind<block_end) - block_size*(data_n-1);
            [spike_ind,spike_clust] = spykingcircus_cleaner(spike_ind,spike_clust);                                
    end

    %%--------------------- Two main values from the detection part: spike_ind, spike_clust + picked_components, picked_comp_top
    %% 3. RAP-MUSIC (2) dipole fitting

    [IndMax_n, ValMax_n, ind_m_n, spikeind_n] = ...
     spike_localization(spike_ind, ...
                        Data, ...
                        G3, ...
                        channel_type, ...
                        parameters.rap_music.f_low_RAP, ...
                        parameters.rap_music.f_high_RAP, ...
                        parameters.rap_music.spikydata, ...
                        picked_components, ...
                        picked_comp_top, ...
                        parameters.corr_thresh, ...
                        parameters.rap_music.RAP, ...
                        '','');
    l = size(IndMax);
    IndMax   = [IndMax IndMax_n +block_size*(data_n-1)];
    ValMax   = [ValMax ValMax_n +block_size*(data_n-1)];
    spikeind = [spikeind; spikeind_n + block_size*(data_n-1)];
    ind_m    = [ind_m  ind_m_n  + l(2)];
end
save(paths.sources_saving_path, 'IndMax','ValMax','ind_m','spikeind')
%% 4. Clustering
clear cluster
if parameters.spikes_detection == 1 % for manual spikes
    %Nmin = 1;
    if size(IndMax) ~= size(spike_ind)
        spike_ind  = spike_ind';
    end
    cluster{1,1} = [IndMax ; spike_ind  ; ValMax ; 1:length(IndMax); zeros(size(spike_ind))];
else
    % cluster creation space
    %try
        cluster = clustering(spike_ind, ... 
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
        
        if parameters.spikes_detection == 3
            % refine clusters throwing away multiple detection of spikes, only for Spyking Circus
            cluster = spykingcircus_cleaner_aftecluster(cluster);
        end
    %catch
    %    if size(IndMax) ~= size(spike_ind)
    %        spike_ind  = spike_ind';
    %    end
    %    cluster{1,1} = [IndMax ; spike_ind  ; ValMax ; 1:length(IndMax)];
    %end
    
end

save(paths.clusters_saving_path ,'cluster')
save(paths.parameters_saving_path ,'parameters')
end

