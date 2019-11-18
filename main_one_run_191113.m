% -------------------------------------------------------------------------
% All steps, one case
% -------------------------------------------------------------------------
% INPUTS:   see main_aspire_all_19113.m
%
% OUTPUTS:

% 1. ICA
%     ICA_spikes_mat = strcat([resultsdir_root, subj_name, '\ASPIRE\',  'ICA_based_', file_name_short, '_' channel_type '.mat']);
%     save(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
% 2. dipole_fit MUSIC
%     RapMusicOutputMAt =  [resultsdir_root, subj_name, '\ASPIRE\sources_'  spikes_extraction '_' channel_type '.mat']
%     save(RapMusicOutputMAt,    'IndMax','ValMax','ind_m','spikeind')
% 3. clustering
%     cluster_out  = cluster_out_spycirc(cluster, G3);
%     csvwrite([resultsdir_root, subj_name,  '\ASPIRE\cluster_out_' spikes_extraction '_' channel_type '.csv'],cluster_out);
% 4. single spike pics
% 5. global pic for each cluster
%     saveas(gcf,[save_param.resultsdir_root, save_param.subj_name, save_param.results_subfolder '\Clusters_' save_param.spikes_extraction '_' channel_type '_' num2str(clust_num)  '.fig'])
%     saveas(gcf,[save_param.resultsdir_root, save_param.subj_name, save_param.results_subfolder '\Clusters_' save_param.spikes_extraction '_' channel_type '_' num2str(clust_num)  '.png'])

    

%
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com


%% 1. Export everything from brainstorm

for channel_type_loop =  2
    switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'mag';
        case 2, channel_type = 'grad';
    end
    labels = labels(3:3:end); % for plotting with fieldtrip

    %% 2. Spike detection
    for spikes_detection = detection_type%1:3
        try
            switch spikes_detection
                case 1 % visual markings
                    spikes_extraction = 'visual';
                    try
                    manual_data = csvread(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                        cases_unique_for_anat, '_', file_name_short, '\Aspire\', ...
                        'Manual_spikes_', file_name_short, '.csv']));
                    picked_components = [];
                    picked_comp_top = [];
                    spike_ind = manual_data(:,1);
                    spike_ind = spike_ind(spike_ind<600-30)*1000;
                    spcirc_clust = [];
                    catch
                       continue  
                    end
                    
                    
                case 2 % ICA based
                    spikes_extraction = 'ICA_based';
                    decision = 0.9; % the amplitude threshold for decision
                    f_low = 3; % bandpass filter before the ICA decomposition
                    f_high = 70;
                    ICA_spikes_mat = strcat([resultsdir_root, subj_name, '\ASPIRE\', ...
                        'ICA_based_', file_name_short, '_' channel_type '.mat']);
                    if newdataset
                        [spike_ind, picked_components, picked_comp_top, component_indicatior] = ...
                            ICA_detection(Data, G3, channel_type, decision, f_low, f_high);
                        save(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
                    end
                    load(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
                    spcirc_clust = [];
                    spike_clust = zeros(size(spike_ind))';
                    spike_clust = spike_clust(spike_ind<600000-30 & spike_ind>41);
                    spike_ind = spike_ind(spike_ind<600000-30 & spike_ind>41);
                    
                    
                case 3 % Spiking circus based
                    spikes_extraction = 'SpyCir_based';
                    spcirc_data = csvread(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                        cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'Templates_',...
                        file_name_short,'_', channel_type, '.csv']),1,0);
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
            if computation_source | newdataset
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
                xlim([0.95 1])
                
                

                RapMusicOutputMAt =  [resultsdir_root, subj_name, '\ASPIRE\sources_'  spikes_extraction '_' channel_type '.mat']
                save(RapMusicOutputMAt,...
                    'IndMax','ValMax','ind_m','spikeind')
                
            end
            load(RapMusicOutputMAt,...
                'IndMax','ValMax','ind_m','spikeind')
            
            fig_gof_name = [spikes_extraction, '_', channel_type];
            figure('Name','fig_gof_name','visible','off')
            histogram(ValMax)
            saveas(gcf,[resultsdir_root, subj_name, ...
                '\ASPIRE\GOF_hist_',spikes_extraction, '_', channel_type, '.bmp']);
            close(gcf);
            
            %% 4. Clustering
            if computation_clusters
                
                % load dipoles
                 load(RapMusicOutputMAt,...
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
                
                % Save events 1 - ALL DETECTED SPIKES
                events = spike_ind/1000;
                save(strcat([resultsdir_root, subj_name, '\ASPIRE\',...
                    'EVENTS1_' spikes_extraction '_' file_name_short, '_', channel_type, '.mat']),'events');
                gof_for_events = ValMax;
                save(strcat([resultsdir_root, subj_name, '\ASPIRE\',...
                    'GOF1_' spikes_extraction '_' file_name_short, '_', channel_type, '.mat']),'gof_for_events');
                
                % Save events 2 - ALL SPIKES DETECTED and FITTED above subcorr THRESHOLD 
                events = spike_ind(ind_m)/1000;
                save(strcat([resultsdir_root, subj_name, '\ASPIRE\',...
                    'EVENTS2_' spikes_extraction '_' file_name_short, '_', channel_type, '.mat']),'events');
                % End saving events
                
                disp(['Subcorr threshold: ', num2str(corr_thresh), ' Number of spike found: ', ...
                    num2str(length(ind_m))]);
                
                thr_dist = THR_DIST; % maximal distance from the center of the cluster (radius) in m
                Nmin     = N_MIN; % minimum number of sources in one cluster
                
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
                try
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
                    save_param.results_subfolder  = ['\ASPIRE\clusters_PICS'];
                    save_param.spikes_extraction  = spikes_extraction;
                    mkdir([resultsdir_root subj_name save_param.results_subfolder])
                    plot_clusters_191015(Data, channel_type, f_low_vis, f_high_vis, cortex, ...
                        spike_trials, maxamp_spike_av, spike_ts, cluster, ...
                        channels, G3, MRI, corr_thresh, save_param ) %epi_plot_autoALLCLUSTERS
                    
                catch
                    draw_and_save_plots_err = 'Error';
                end
            end
            
            
            %% 7. all to all cluster propagation probability
            if 0
                time_w  = .1; % window of interest is  [-time_w time_w] in seconds
                distr = all2all_prop(cluster,time_w)
            end
            
            
            %% 8. write clusters in csv to file
            if computation_clusters
                cluster_out  = cluster_out_spycirc(cluster, G3);
                csvwrite([resultsdir_root, subj_name,  '\ASPIRE\cluster_out_' spikes_extraction '_' channel_type '.csv'],cluster_out);
                
                events = cluster_out(:,1)/1000+offset_time;
                save(strcat([resultsdir_root, subj_name,  '\ASPIRE\',...
                    'EVENTS3_' spikes_extraction '_' file_name_short, '_', channel_type, '.mat']),'events');
                
            end
            
            %% plot spikes of same cluster
            for cl = 1:length(cluster)
                
                spikes_fitted = [cluster_out(:,1) cluster_out(:,1) cluster_out(:,3) cluster_out(:,7)  cluster_out(:,2)]%  [record time in s for brainstorm ¦ idx detecton sample  ¦ subspace corr ¦ leadfield_orig_index ¦ cluster]
                spikes_fitted = sort(spikes_fitted,1)
                spikes_fitted(:,1) = spikes_fitted(:,1)/1000+offset_time % cluster_out_times_only
                
                spikes_fitted = spikes_fitted(find(spikes_fitted(:,end) == cl),:);
                path_spikes_plot = [resultsdir_root, subj_name '\ASPIRE\SPIKES\cluster' num2str(cl) '\']
                mkdir(path_spikes_plot)
                TOPO = plot_spikes_ER_TOPO(spikes_fitted, path_spikes_plot, Data, channel_type, f_low, f_high,  G3,  channels)
                
                %write cluster by cluster to bs event file
                events = spikes_fitted(:,1);
                save ([resultsdir_root, subj_name '\ASPIRE\EVENTS3_cluster' num2str(cl) ],'events')
                
                
            end
            
            %save all fitted spikes
            spikes_fitted = [cluster_out(:,1) cluster_out(:,1) cluster_out(:,3) cluster_out(:,7)  cluster_out(:,2)]%  [record time in s for brainstorm ¦ idx detecton sample  ¦ subspace corr ¦ leadfield_orig_index ¦ cluster]
            spikes_fitted = sort(spikes_fitted,1)
            spikes_fitted(:,1) = spikes_fitted(:,1)/1000+offset_time % cluster_out_times_only
            save ([resultsdir_root, subj_name '\ASPIRE\spikes_fitted' ],'spikes_fitted')
            
            %% saving
            if draw_and_save_plots2
                try
                    f_low_vis  = 2; % bandpass filter for visualization
                    f_high_vis = 50;
                    %% scatter
                    figure('Name','Clusters raster plot','visible','on')
                    scatter(cluster_out(:,1)/1000,cluster_out(:,2),'o')
                    grid minor
                    
                    param.subj_name             = subj_name;
                    param.spikes_detection   = spikes_detection;
                    param.spikes_extraction  = spikes_extraction;
                    param.channel_type         = channel_type;
                    param.Nmin                       = Nmin;
                    param.thr_dist                   = thr_dist;
                    param.corr_thresh            = corr_thresh;
                    param.f_low_RAP            = f_low_RAP;
                    param.f_high_RAP           = f_high_RAP;
                    param.f_low_vis                = f_low_vis; % bandpass filter for visualization
                    param.f_high_vis              = f_high_vis;
                    param.time_w                   = time_w;
                    param.distr                        = distr;
                    
                    save([resultsdir_root, subj_name '\ASPIRE\results_'  spikes_extraction '_' channel_type '.mat'],...
                        'cluster','param')
                    
                    % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder '\Clusters_' spikes_extraction '_' channel_type '.bmp'])
                    % saveas(fig_cluster,[resultsdir_root, subj_name, results_subfolder '\Clusters_' spikes_extraction '_' channel_type '.fig'])
                catch
                    draw_and_save_plots2_err = 'Error';
                end
            end
            
        catch
            source_clusters_err   = 'Error';
        end
    end
    
    %% Plot all manual spikes -----------0
    if plot_all_manual_spikes
        try
            spike_time_TP = load([resultsdir_root, subj_name, results_subfolder,...
                'cluster_out_visual_', channel_type '.csv']);
            
            
            spike_time_ICA_ev1 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS1_',...
                'ICA_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_ICA_gof = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'GOF1_',...
                'ICA_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_ICA  = [spike_time_ICA_ev1.events'*1000 spike_time_ICA_gof.gof_for_events'];
            spike_time_ICA_ev2 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS2_',...
                'ICA_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_ICA_ev2 = spike_time_ICA_ev2.events'*1000;
            spike_time_ICA_ev3 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS3_',...
                'ICA_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_ICA_ev3 = spike_time_ICA_ev3.events*1000;
            
            spike_time_SPC_ev1 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS1_',...
                'SpyCir_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_SPC_gof = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'GOF1_',...
                'SpyCir_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_SPC  = [spike_time_SPC_ev1.events*1000 spike_time_SPC_gof.gof_for_events'];
            spike_time_SPC_ev2 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS2_',...
                'SpyCir_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_SPC_ev2 = spike_time_SPC_ev2.events*1000;
            spike_time_SPC_ev3 = load(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...
                cases_unique_for_anat, '_', file_name_short, '\Aspire\', 'EVENTS3_',...
                'SpyCir_based_', file_name_short, '_', channel_type ,'.mat']));
            spike_time_SPC_ev3 = spike_time_SPC_ev3.events*1000;
            
            f_low_vis  = 2;
            f_high_vis = 50;
            
            mkdir([resultsdir_root subj_name '\' results_subfolder 'Manual_spikes_FN_grad'])
            mkdir([resultsdir_root subj_name '\' results_subfolder 'Manual_spikes_TP_grad'])
            mkdir([resultsdir_root subj_name '\' results_subfolder 'Manual_spikes_FN_mag'])
            mkdir([resultsdir_root subj_name '\' results_subfolder 'Manual_spikes_TP_mag'])
            path = [resultsdir_root subj_name '\' results_subfolder 'Manual_spikes_'];
            
            
            plot_manual_spikes_191015(spike_time_TP, path, Data, f_low_vis, f_high_vis, ...
                channel_type, G3, channels, spike_time_ICA, spike_time_SPC, ...
                spike_time_ICA_ev2,  spike_time_ICA_ev3, spike_time_SPC_ev2, ...
                spike_time_SPC_ev3)
            
        catch
            plot_all_manual_spikes_err = 'Error';
        end
    end
    
    
    %% Errors list
    xlsTAB = subj_name;
    backsalsh = find(xlsTAB=='\');
    if backsalsh
        xlsTAB(backsalsh)='_';
    end
    
    fname = [resultsdir_root 'Aspire_ROC\Errors_COR_TR_' num2str(CORR_THR) '.xlsx'];
    xlswrite(fname,{'MAG'},xlsTAB,'B1');
    xlswrite(fname,{'GRAD'},xlsTAB,'C1');
    xlswrite(fname,{'source_clusters'},xlsTAB,'A2');
    xlswrite(fname,{'draw_and_save_plots'},xlsTAB,'A3');
    xlswrite(fname,{'draw_and_save_plots2'},xlsTAB,'A4');
    xlswrite(fname,{'plot_all_manual_spikes'},xlsTAB,'A5');
    
    if strcmp(channel_type, 'mag') == 1
        xlswrite(fname,{source_clusters_err},xlsTAB,'B2');
        xlswrite(fname,{draw_and_save_plots_err},xlsTAB,'B3');
        xlswrite(fname,{draw_and_save_plots2_err},xlsTAB,'B4');
        xlswrite(fname,{plot_all_manual_spikes_err},xlsTAB,'B5');
    else
        xlswrite(fname,{source_clusters_err},xlsTAB,'C2');
        xlswrite(fname,{draw_and_save_plots_err},xlsTAB,'C3');
        xlswrite(fname,{draw_and_save_plots2_err},xlsTAB,'C4');
        xlswrite(fname,{plot_all_manual_spikes_err},xlsTAB,'C5');
    end
    source_clusters_err   = '';
    draw_and_save_plots_err = '';
    draw_and_save_plots2_err = '';
    plot_all_manual_spikes_err = '';
    
end

%% Plot BIGPIC
if plot_big_pic
    try
        BIGPIC_clusters_191016(subj_name, results_subfolder, ...
            resultsdir_root, cortex, CORR_THR)
        
    catch
        plot_big_pic_err = 'Error';
    end
end


%% ROC curves
if computation_ROC
    try
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
    catch
        computation_ROC_err = 'Error';
    end
end

%% Errors list

fname = [resultsdir_root 'Aspire_ROC\Errors_COR_TR_' num2str(CORR_THR) '.xlsx'];

xlswrite(fname,{'computation_ROC'},xlsTAB,'A6');
xlswrite(fname,{computation_ROC_err},xlsTAB,'D6');
xlswrite(fname,{'plot_big_pic'},xlsTAB,'A9');
xlswrite(fname,{plot_big_pic_err},xlsTAB,'D9');

computation_ROC_err = '';
plot_big_pic_err = '';



