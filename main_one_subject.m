function main_one_subject(cortex, Data, G3, paths, parameters)

% -------------------------------------------------------------------------
% All steps, one case
% -------------------------------------------------------------------------
% INPUTS:
% detection_type -- 1:visual, 2:ICA-based, 3:Spyking Circus
%
% PATHS
%
% path_vis_detections -- path to the visual detections file (csv)
% path_ICA_detections -- path to the ICA detections (mat)
% path_SPC_detections -- path to the Spyking Circus detections (csv)
% resultsdir_root
% subj_name
% results_subfolder
% 
% PARAMETERS
%
% mute_mode -- not show plots
% computation_source
% computation_clusters
% draw_and_save_plots
% draw_and_save_plots2
% computation_ROC
% plot_big_pic
% Data
% G3
% THR_DIST - maximal distance from the center of the cluster (radius) in m
% N_MIN - minimum number of sources in one cluster
% roc_xlsx_fname - path to the main ROC saving file
% roc_labels_xlsx_fname - path to the labels ROC saving file
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
                picked_components = [];
                picked_comp_top = [];
                spike_ind = manual_data(:,1);
                spike_ind = spike_ind(spike_ind<600-30)*1000;
                spcirc_clust = [];
                
                
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
                spcirc_clust = [];
                spike_clust = zeros(size(spike_ind))';
                spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                
                
            case 3 % Spiking circus based
                spikes_extraction = parameters.detection.SPC.spikes_extraction;
                spcirc_data = csvread([paths.path_SPC_detections,'_', channel_type, '.csv'],1,0);
                picked_components = [];
                picked_comp_top = [];
                spike_ind = spcirc_data(:,1);
                spike_clust = spcirc_data(:,3);
                spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                [spike_ind,spike_clust] = spykingcircus_cleaner(spike_ind,spike_clust);
                
                
        end
        close all
        
        %% 3. RAP-MUSIC (2) dipole fitting
        if parameters.computation_source            
            if spikes_detection ~= 1
                corr_thresh = 0.0;
            else
                corr_thresh = parameters.CORR_THR; %0.95
            end
            
            % corr_thresh = back to quantile(ValMax, 0.95)
            [IndMax, ValMax, ~, spikeind] = spike_localization(spike_ind, Data, G3, ...
                channel_type, parameters.rap_music.f_low_RAP, parameters.rap_music.f_high_RAP, ...
                parameters.rap_music.spikydata, picked_components, ...
                picked_comp_top, corr_thresh, parameters.rap_music.RAP);
            
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

            % set  CORR_TRESH
            if spikes_detection ~= 1
                corr_thresh = prctile(ValMax,85);
            else
                corr_thresh = parameters.CORR_THR; %0.95
                %corr_thresh = CORR_THR;
            end
            ind_m = find((ValMax > corr_thresh));
            
            %disp(['Subcorr threshold: ', num2str(corr_thresh), ' Number of spike found: ', ...
            %    num2str(length(ind_m))]);
            
            clear cluster
            if spikes_detection == 1 % for manual spikes
                %Nmin = 1;
                if size(IndMax) ~= size(spike_ind)
                    spike_ind  = spike_ind';
                end
                cluster{1,1} = [IndMax ; spike_ind  ; ValMax ; 1:length(IndMax)];
            else
                % cluster creation space
                try
                    cluster = clustering(spike_ind, G3, parameters.clustering.N_MIN, ValMax, IndMax, ind_m, ...
                        parameters.clustering.THR_DIST, 1, cortex, parameters.rap_music.RAP, spikeind, spike_clust);
                    
                    if spikes_detection == 3
                        % refine clusters throwing away multiple detection of spikes
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
            cluster_out  = cluster_out(cluster, G3);
            
            csvwrite([paths.path_cluster_out spikes_extraction '_' channel_type '.csv'], cluster_out);
            
        end
        
        %% 5. Activation on sources
        if parameters.draw_and_save_plots
            
            [spike_trials, maxamp_spike_av, ~, spike_ts] = ...
                source_reconstruction(Data, G3, channel_type, cluster, ...
                parameters.draw.f_low, parameters.draw.f_high);
            
            close all
            %% 6. Big plot
            
            save_param.resultsdir_root    = paths.resultsdir_root;
            save_param.subj_name          = paths.subj_name;
            save_param.results_subfolder  = [paths.results_subfolder '\clusters'];
            save_param.spikes_extraction  = spikes_extraction;
            
            plot_clusters_191015(Data, channel_type, parameters.draw.f_low_vis, ...
                parameters.draw.f_high_vis, cortex, ...
                spike_trials, maxamp_spike_av, spike_ts, cluster, ...
                channels, G3, MRI, corr_thresh, save_param ) %epi_plot_autoALLCLUSTERS
        end
        
        
        %% 7. all to all cluster propagation probability
        if 0
            time_w  = .1; % window of interest is  [-time_w time_w] in seconds
            distr = all2all_prop(cluster,time_w)
        end
        
        %% saving
        if parameters.draw_and_save_plots2
            %% scatter
            figure('Name','Clusters raster plot','visible','off')
            scatter(cluster_out(:,1)/1000,cluster_out(:,2),'.')
            grid minor
            
            param.subj_name             = subj_name;
            param.spikes_detection      = spikes_detection;
            param.spikes_extraction     = spikes_extraction;
            param.channel_type          = channel_type;
            param.Nmin                  = parameters.clustering.N_MIN;
            param.thr_dist              = parameters.clustering.THR_DIST;
            param.corr_thresh           = corr_thresh;
            %             param.f_low_RAP             = f_low_RAP;
            %             param.f_high_RAP            = f_high_RAP;
            param.f_low_vis             = parameters.draw.f_low_vis; % bandpass filter for visualization
            param.f_high_vis            = parameters.draw.f_high_vis;
            %             param.time_w                = time_w;
            %             param.distr                 = distr;
            
            save([results_saving_path spikes_extraction '_' channel_type '.mat'] ,'cluster','param')
            
            % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder, ...
            %        '\Clusters_' spikes_extraction '_' channel_type '.bmp'])
            % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder, ...
            %        '\Clusters_' spikes_extraction '_' channel_type '.fig'])
            
        end
        
    end
    
    %% Plot BIGPIC
    if parameters.plot_big_pic
        plot_bigpic(paths.subj_name, paths.results_saving_path, cortex, paths.bigpic_saving_path)
        
    end
    
    %% ROC curves
    if parameters.computation_ROC
        
        % Load Spyking Circus detected timestamps
        SPC_grad = load([paths.path_cluster_out 'SpyCir_based_grad.csv']);
        SPC_mag = load([paths.path_cluster_out 'SpyCir_based_mag.csv']);
        
        % Load ICA detected timestamps
        ICA_grad = load([paths.path_cluster_out 'ICA_based_grad.csv']);
        ICA_mag = load([paths.path_cluster_out 'ICA_based_mag.csv']);
        
        % Load visual timestamps
        visual = load([paths.path_cluster_out 'visual_grad.csv']);
        
        % Compute and save all results in the excel file
        ROC(parameters.detection_type, ICA_grad, ICA_mag, SPC_grad, SPC_mag, visual, ...
            cortex, paths.roc_xlsx_fname, paths.roc_labels_xlsx_fname)
        
    end
    
    
end


