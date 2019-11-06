function = main_one_subject(detection_type, path_vis_detections, path_ICA_detections, path_SPC_detections,
newdataset, resultsdir_root, subj_name, results_subfolder, mute_mode, computation_source, computation_clusters,
draw_and_save_plots, raw_and_save_plots2, CORR_THR, Data, G3)

% -------------------------------------------------------------------------
% All steps, one case
% -------------------------------------------------------------------------
% INPUTS:
% detection_type -- 1:visual, 2:ICA-based, 3:Spyking Circus
% path_vis_detections -- path to the visual detections file (csv)
% path_ICA_detections -- path to the ICA detections (mat)
% path_SPC_detections -- path to the Spyking Circus detections (csv)
% resultsdir_root
% subj_name
% results_subfolder
% mute_mode -- not show plots
% computation_source
% computation_clusters
% draw_and_save_plots
% raw_and_save_plots2
% Data
% G3
% computation_ROC
% plot_big_pic
%
%
% OUTPUTS:
%
% _______________________________________________________
%
for channel_type_loop = 1:2
    switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'mag';
        case 2, channel_type = 'grad';
    end
    
    %% 2. Spike detection
    for spikes_detection = detection_type
        
        switch spikes_detection
            case 1 % visual markings
                spikes_extraction = 'visual';
                manual_data = csvread(path_vis_detections);
                picked_components = [];
                picked_comp_top = [];
                spike_ind = manual_data(:,1);
                spike_ind = spike_ind(spike_ind<600-30)*1000;
                spcirc_clust = [];
                
                
            case 2 % ICA based
                spikes_extraction = 'ICA_based';
                decision = 0.9; % the amplitude threshold for decision
                f_low = 3; % bandpass filter before the ICA decomposition
                f_high = 70;
                ICA_spikes_mat = strcat([path_ICA_detections,'_', channel_type, '.mat']);
                
                if newdataset
                    [spike_ind, picked_components, picked_comp_top] = ...
                        ICA_detection(Data, G3, channel_type, decision, f_low, f_high);
                    save(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top')
                end
                
                load(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top')
                spcirc_clust = [];
                spike_clust = zeros(size(spike_ind))';
                spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                
                
            case 3 % Spiking circus based
                spikes_extraction = 'SpyCir_based';
                spcirc_data = csvread(strcat([path_SPC_detections,'_', channel_type, '.csv'])),1,0);
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
        if computation_source
            f_low_RAP  = 10;
            f_high_RAP = 200;
            spikydata = 0;
            %RAP = 'RAP'; corr_thresh = 0.99;
            RAP = 'not';
            if spikes_detection ~= 1
                corr_thresh = 0.0;
            else
                corr_thresh = CORR_THR; %0.95
            end
            % corr_thresh = back to quantile(ValMax, 0.95)
            [IndMax, ValMax, ind_m, spikeind] = spike_localization(spike_ind, Data, G3, ...
                channel_type, f_low_RAP, f_high_RAP, spikydata, picked_components, ...
                picked_comp_top, corr_thresh, RAP);
            
            save([resultsdir_root, subj_name, results_subfolder '\sources_'  spikes_extraction '_' channel_type '.mat'],...
                'IndMax','ValMax','ind_m','spikeind')
        end
        
        load([resultsdir_root, subj_name, results_subfolder '\sources_'  spikes_extraction '_' channel_type '.mat'],...
            'IndMax','ValMax','ind_m','spikeind')
        
        fig_gof_name = [spikes_extraction, '_', channel_type];
        
        figure('Name','fig_gof_name','visible','off')
        histogram(ValMax)
        saveas(gcf,[resultsdir_root, subj_name, results_subfolder,...
            'GOF_hist_',spikes_extraction, '_', channel_type, '.bmp']);
        close(gcf);
        
        %% 4. Clustering
        if computation_clusters
            
            % load dipoles
            load([resultsdir_root, subj_name, results_subfolder '\sources_'  spikes_extraction '_' channel_type '.mat'],...
                'IndMax','ValMax','ind_m','spikeind')
            RAP = 'not';
            % set  CORR_TRESH
            if spikes_detection ~= 1
                corr_thresh = prctile(ValMax,85);
            else
                corr_thresh = CORR_THR; %0.95
                %corr_thresh = CORR_THR;
            end
            ind_m = find((ValMax > corr_thresh));
            
            
            
            disp(['Subcorr threshold: ', num2str(corr_thresh), ' Number of spike found: ', ...
                num2str(length(ind_m))]);
            
            thr_dist = THR_DIST; % maximal distance from the center of the cluster (radius) in m
            Nmin = N_MIN; % minimum number of sources in one cluster
            
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
                    cluster = clustering(spike_ind, G3, Nmin, ValMax, IndMax, ind_m, ...
                        thr_dist, 1,cortex, RAP, spikeind, spike_clust);
                    
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
                % cluster creation labels
                %             [~, cluster] = clustering_by_labels(spike_ind, ValMax, IndMax,...
                %                 ind_m,RAP, spikeind, spike_clust, cortex.Atlas);
                
                
            end
            close all
        end
        
        %% 5. Activation on sources
        if draw_and_save_plots
            
            f_low  = 3;
            f_high = 50;
            [spike_trials, maxamp_spike_av, channels_maxamp, spike_ts] = ...
                source_reconstruction(Data, G3, channel_type, cluster, ...
                f_low, f_high);
            
            close all
            %% 6. Big plot
            f_low_vis  = 2; % bandpass filter for visualization
            f_high_vis = 50;
            
            %         epi_plot(Data, channel_type, f_low_vis, f_high_vis, cortex, ...
            %             spike_trials, maxamp_spike_av, spike_ts, cluster, ...
            %             channels, G3, MRI, corr_thresh)
            
            save_param.resultsdir_root    =  resultsdir_root;
            save_param.subj_name          = subj_name;
            save_param.results_subfolder  = [results_subfolder '\clusters'];
            save_param.spikes_extraction  = spikes_extraction;
            
            plot_clusters_191015(Data, channel_type, f_low_vis, f_high_vis, cortex, ...
                spike_trials, maxamp_spike_av, spike_ts, cluster, ...
                channels, G3, MRI, corr_thresh, save_param ) %epi_plot_autoALLCLUSTERS
            
            
        end
        
        
        %% 7. all to all cluster propagation probability
        if 0
            time_w  = .1; % window of interest is  [-time_w time_w] in seconds
            distr = all2all_prop(cluster,time_w)
        end
        
        
        %% 8. write clusters in csv to file
        if computation_clusters
            cluster_out  = cluster_out_spycirc(cluster, G3);
            csvwrite([resultsdir_root, subj_name, results_subfolder '\cluster_out_' spikes_extraction '_' channel_type '.csv'],cluster_out);
            
        end
        
        %% saving
        if draw_and_save_plots2
            
            f_low_vis  = 2; % bandpass filter for visualization
            f_high_vis = 50;
            %% scatter
            figure('Name','Clusters raster plot','visible','off')
            scatter(cluster_out(:,1)/1000,cluster_out(:,2),'.')
            grid minor
            
            param.subj_name             = subj_name;
            param.spikes_detection      = spikes_detection;
            param.spikes_extraction     = spikes_extraction;
            param.channel_type          = channel_type;
            param.Nmin                  = Nmin;
            param.thr_dist              = thr_dist;
            param.corr_thresh           = corr_thresh;
            %             param.f_low_RAP             = f_low_RAP;
            %             param.f_high_RAP            = f_high_RAP;
            param.f_low_vis             = f_low_vis; % bandpass filter for visualization
            param.f_high_vis            = f_high_vis;
            %             param.time_w                = time_w;
            %             param.distr                 = distr;
            
            save([resultsdir_root, subj_name, results_subfolder '\results_'  spikes_extraction '_' channel_type '.mat'],...
                'cluster','param')
            
            % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder '\Clusters_' spikes_extraction '_' channel_type '.bmp'])
            % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder '\Clusters_' spikes_extraction '_' channel_type '.fig'])
            
        end
        
    end
end

%% Plot BIGPIC
if plot_big_pic
    BIGPIC_clusters_191016(subj_name, results_subfolder, ...
        resultsdir_root, cortex, CORR_THR)
    
end
%% ROC curves
if computation_ROC
    
    xlsTAB = subj_name;
    backsalsh = find(xlsTAB=='\');
    if backsalsh
        xlsTAB(backsalsh)='_';
    end
    
    for spikes_detection = detection_type%1:3
        
        switch spikes_detection
            
            case 2
                
                % ICA
                fname = [resultsdir_root 'Aspire_ROC\ROC_COR_TR_' num2str(CORR_THR) '.xlsx'];
                fname_labels = [resultsdir_root 'Aspire_ROC\Labels_COR_TR_' num2str(CORR_THR) '.xlsx'];
                
                ICA_grad = load([resultsdir_root, subj_name, results_subfolder,...
                    'cluster_out_ICA_based_grad.csv']);
                ICA_mag = load([resultsdir_root, subj_name, results_subfolder,...
                    'cluster_out_ICA_based_mag.csv']);
                ICA_visual = load([resultsdir_root, subj_name, results_subfolder,...
                    'cluster_out_visual_grad.csv']);
                ICA_visual_mag = load([resultsdir_root, subj_name, results_subfolder,...
                    'cluster_out_visual_grad.csv']);
                ICA_visual_grad = load([resultsdir_root, subj_name, results_subfolder, ...
                    'cluster_out_visual_grad.csv']);
                
                [ICA_labels_results, ICA_roc_labels] = roc_curve_labels(ICA_visual_mag, ICA_visual_grad, ICA_mag,  ICA_grad, cortex.Atlas);
                ICA_roc = roc_curve(ICA_visual(ICA_visual(:,1)<600000,:), ICA_grad, ICA_mag);
                ICA_roc_spatial = roc_curve_spatial(ICA_visual_grad, ICA_visual_mag, ICA_grad, ICA_mag);
                
                xlswrite(fname,ICA_roc,xlsTAB,'A6');
                xlswrite(fname,{'ICA'},xlsTAB,'A6');
                xlswrite(fname,ICA_roc_spatial,xlsTAB,'A16');
                xlswrite(fname,{'ICA'},xlsTAB,'A16');
                xlswrite(fname,ICA_roc_labels,xlsTAB,'A26');
                xlswrite(fname,{'ICA'},xlsTAB,'A26');
                xlswrite(fname_labels,ICA_labels_results,xlsTAB,'I1');
                xlswrite(fname_labels,{'ICA'},xlsTAB,'I1');
                
            case 3
                % Spyking Circus
                fname = [resultsdir_root 'Aspire_ROC\ROC_COR_TR_' num2str(CORR_THR) '.xlsx'];
                fname_labels = [resultsdir_root 'Aspire_ROC\Labels_COR_TR_' num2str(CORR_THR) '.xlsx'];
                
                SPC_grad = load([resultsdir_root, subj_name,results_subfolder, ...
                    'cluster_out_SpyCir_based_grad.csv']);
                SPC_mag = load([resultsdir_root, subj_name, results_subfolder, ...
                    'cluster_out_SpyCir_based_mag.csv']);
                SPC_visual = load([resultsdir_root, subj_name, results_subfolder, ...
                    'cluster_out_visual_grad.csv']);
                SPC_visual_mag = load([resultsdir_root, subj_name, results_subfolder, ...
                    'cluster_out_visual_grad.csv']);
                SPC_visual_grad = load([resultsdir_root, subj_name, results_subfolder, ...
                    'cluster_out_visual_grad.csv']);
                
                [SPC_labels_results, SPC_roc_labels] = roc_curve_labels(SPC_visual_mag, SPC_visual_grad, SPC_mag,  SPC_grad, cortex.Atlas);
                SPC_roc = roc_curve(SPC_visual(SPC_visual(:,1)<600000,:), SPC_grad, SPC_mag);
                SPC_roc_spatial = roc_curve_spatial(SPC_visual_grad, SPC_visual_mag, SPC_grad, SPC_mag);
                
                xlswrite(fname,SPC_roc,xlsTAB,'A1');
                xlswrite(fname,{'SPC'},xlsTAB,'A1');
                xlswrite(fname,SPC_roc_spatial,xlsTAB,'A11');
                xlswrite(fname,{'SPC'},xlsTAB,'A11');
                xlswrite(fname,SPC_roc_labels,xlsTAB,'A21');
                xlswrite(fname,{'SPC'},xlsTAB,'A21');
                xlswrite(fname_labels,SPC_labels_results,xlsTAB,'A1');
                xlswrite(fname_labels,{'SPC'},xlsTAB,'A1');
                
        end
    end
    
end