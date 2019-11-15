function main_one_subject(cortex, Data, G3, MRI, channels, paths, parameters)

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


for channel_type_loop = 1:2
    switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'mag';
        case 2, channel_type = 'grad';
    end
    
    %% 2. Spike detection
    for spikes_detection = parameters.detection_type
        
        switch spikes_detection
            case 1 % visual markings
                spikes_extraction =  parameters.detection.visual.spikes_extraction;
                manual_data = csvread(paths.path_vis_detections);
                picked_components = []; % ICA relevant
                picked_comp_top = []; % ICA relevant
                spcirc_clust = []; % SPC relevant (maybe delete)
                spike_ind = manual_data(:,1);
                spike_ind = spike_ind(spike_ind<600-30)*1000;
                
            case 2 % ICA based
                spikes_extraction = parameters.detection.ICA.spikes_extraction;
                ICA_spikes_mat_saving_path = [paths.path_ICA_detections,'_', channel_type, '.mat'];
                
                if parameters.newdataset
                    [spike_ind, picked_components, picked_comp_top] = ...
                        ICA_detection(Data, G3, channel_type, parameters.detection.ICA.decision, ...
                        parameters.detection.ICA.f_low, parameters.detection.ICA.f_high);
                    save(ICA_spikes_mat_saving_path,'spike_ind', 'picked_components', 'picked_comp_top')
                end
                
                load(ICA_spikes_mat_saving_path,'spike_ind', 'picked_components', 'picked_comp_top')
                spcirc_clust = []; % SPC relevant (maybe delete)
                spike_clust = zeros(size(spike_ind))';
                spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                
                
            case 3 % Spiking circus based
                spikes_extraction = parameters.detection.SPC.spikes_extraction;
                spcirc_data = csvread([paths.path_SPC_detections,'_', channel_type, '.csv'],1,0);
                picked_components = []; % ICA relevant
                picked_comp_top = []; % ICA relevant
                spike_ind = spcirc_data(:,1);
                spike_clust = spcirc_data(:,3);
                spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                [spike_ind,spike_clust] = spykingcircus_cleaner(spike_ind,spike_clust);
                
                
        end
        close all
        
        %% 3. RAP-MUSIC (2) dipole fitting
        if parameters.computation_source
            
            % corr_thresh = prctile(ValMax,85);
            disp(['Spikes extraction: ' spikes_extraction '; Channels: ' channel_type]);
            
            [IndMax, ValMax, ind_m, spikeind] = spike_localization(spike_ind, Data, G3, ...
                channel_type, parameters.rap_music.f_low_RAP, parameters.rap_music.f_high_RAP, ...
                parameters.rap_music.spikydata, picked_components, picked_comp_top, ...
                spikes_detection, parameters.prctile, parameters.corr_thresh, parameters.rap_music.RAP);
            
            save([paths.sources_saving_path spikes_extraction '_' channel_type '.mat'], ...
                'IndMax','ValMax','ind_m','spikeind')
        else
            load([paths.sources_saving_path spikes_extraction '_' channel_type '.mat'], ...
                'IndMax','ValMax','ind_m','spikeind')
        end
        % Plot and save ValMax
        %        fig_gof_name = [spikes_extraction, '_', channel_type];
        %        figure('Name','fig_gof_name','visible','off')
        %        histogram(ValMax)
        %        saveas(gcf,[resultsdir_root, subj_name, results_subfolder,...
        %            'GOF_hist_',spikes_extraction, '_', channel_type, '.bmp']);
        %        close(gcf);
        
        %% 4. Clustering
        if parameters.computation_clusters
            % load dipoles
            load([paths.sources_saving_path spikes_extraction '_' channel_type '.mat'], ...
                'IndMax','ValMax','ind_m','spikeind')
            
            clear cluster
            if spikes_detection == 1 % for manual spikes
                %Nmin = 1;
                if size(IndMax) ~= size(spike_ind)
                    spike_ind  = spike_ind';
                end
                cluster{1,1} = [IndMax ; spike_ind  ; ValMax ; 1:length(IndMax); zeros(size(spike_ind))];
            else
                % cluster creation space
                try
                    cluster = clustering(spike_ind, G3, parameters.clustering.N_MIN, ValMax, IndMax, ind_m, ...
                        parameters.clustering.THR_DIST, 1, cortex, parameters.rap_music.RAP, spikeind, spike_clust);
                    
                    if spikes_detection == 3
                        % refine clusters throwing away multiple detection of spikes, only for Spyking Circus
                        cluster = spykingcircus_cleaner_aftecluster(cluster);
                    end
                catch
                    if size(IndMax) ~= size(spike_ind)
                        spike_ind  = spike_ind';
                    end
                    cluster{1,1} = [IndMax ; spike_ind  ; ValMax ; 1:length(IndMax)];
                end
                
            end
            
            % Write clusters in csv file
            cluster_out_results = cluster_out(cluster, G3);
            csvwrite([paths.path_cluster_out spikes_extraction '_' channel_type '.csv'], ...
                cluster_out_results);
            
            % save only timestamps
            cluster_out_time_only = cluster_out_results(:,1); % !!!should be in seconds, add first sample, Data.Time(1)
            save([paths.path_cluster_out 'time_only_' spikes_extraction '_' channel_type '.mat'], ...
                'cluster_out_time_only');
        end
        
        %% 5. Activation on sources
        if parameters.draw_and_save_plots && (spikes_detection>1)
            
            [spike_trials, maxamp_spike_av, ~, spike_ts] = ...
                source_reconstruction(Data, G3, channel_type, cluster, ...
                parameters.draw.f_low, parameters.draw.f_high);
            
            close all
            %% 6. Big plot
            
            plot_clusters(Data, channel_type, spikes_extraction, parameters.draw.f_low_vis, ...
                parameters.draw.f_high_vis, cortex, ...
                spike_trials, maxamp_spike_av, spike_ts, cluster, ...
                channels, G3, MRI, prctile(ValMax,parameters.prctile), ...
                parameters.draw.save_clusters, paths.save_cluster_plots, ...
                parameters.mute_mode) %epi_plot_autoALLCLUSTERS
        end
        
        
        %% 7. all to all cluster propagation probability
        if parameters.propagation_probability
            time_w  = .1; % window of interest is  [-time_w time_w] in seconds
            distr = all2all_prop(cluster,time_w);
        end
        
        %% saving results and parameters
        if parameters.save_results
            save([paths.results_saving_path spikes_extraction '_' channel_type '.mat'],'cluster','parameters')
        end
        
    end
end

%% Plot BIGPIC
if parameters.plot_big_pic
    
    plot_bigpic(paths.subj_name, paths.results_saving_path, cortex, paths.bigpic_saving_path)
    
end

%% ROC curves
if parameters.computation_ROC
    
    % Compute and save all results in the excel file
    ROC(paths.subj_name, parameters.detection_type, paths.path_cluster_out, ...
        cortex, paths.roc_xlsx_fname, paths.roc_labels_xlsx_fname)
    
end

%% Overlap between detections
if parameters.compute_overlap
    ICA_and_SPC_spikes(paths.path_cluster_out)
 
    aspire_and_spc_clusters([paths.path_cluster_out 'overlap_mag.mat'], ...
        paths.overlap_saving_path, 'mag', parameters.draw.f_low_vis, ...
        parameters.draw.f_high_vis, Data, channels, 3)
    aspire_and_spc_clusters([paths.path_cluster_out 'overlap_grad.mat'], ...
        paths.overlap_saving_path, 'grad', parameters.draw.f_low_vis, ...
        parameters.draw.f_high_vis, Data, channels, 5)
end
end

