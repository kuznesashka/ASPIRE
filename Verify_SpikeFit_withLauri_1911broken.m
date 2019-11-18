%% all necessary code to start

%% libraries
START_ASPIRE

%% 0. PAREMETERS

%load cases names
load cases_files_191022.mat 
load param_matrix190809.mat

% OUTPUTS
OUTPUT_visclust = []; % how many spikes survive after clustering
ROC_line_excel_spatial_ICA = [ ];
ROC_line_excel_labels_ICA  = [ ];
ROC_line_excel_spatial_SPC = [ ];
ROC_line_excel_labels_SPC  = [ ];
ROC_line_excel_time_ICA    = [ ];
ROC_line_excel_time_SPC    = [ ];

%% Parameters iteration
case_n = 1%[1:4 7:length(cases_files.file_names_short)] %1:4

newdataset = 0;
%
%         for case_n = [11:length(file_names_short)        ]

STAT_ICA_temp = [];
STAT_ICA_spac = [];
STAT_ICA_labe  = [];

STAT_SPC_temp = [];
STAT_SPC_spac = [];
STAT_SPC_labe  = [];

c_tr = [0.955] %var1 = 19 %:length(param_matrix)

var1 = 19;
CORR_THR = c_tr; %param_matrix(var1,1);
THR_DIST = param_matrix(var1,2);
N_MIN    = param_matrix(var1,3);

close all
%         try
subj_name = cases_files.cases{case_n};
cases_unique_for_anat = cases_files.cases_unique{case_n};
%protocol_dir = '\\172.16.118.134\TFedele\Tommaso\Matlab Toolboxes\EPILEPSY\MEG_Tommaso\';
protocol_dir = [hdisk 'Valerii\EPILEPSY\MEG_Tommaso\'];
file_name = cases_files.file_names{case_n};
file_name_short = cases_files.file_names_short{case_n};
resultsdir_root = [hdisk 'Valerii\45_cases\'];
results_subfolder = ['\ASPIRE\CORR_THR_' num2str(CORR_THR) '\'];
% mkdir([resultsdir_root subj_name results_subfolder])
% mkdir([resultsdir_root subj_name results_subfolder 'clusters\'])

computation_source       = 1; % compute dipoles
computation_clusters     = 1; % compute clustering
draw_and_save_plots    = 0; % plot clusters
draw_and_save_plots2  = 1; % save clustering
plot_all_manual_spikes  = 0; % plot all manual spikes
plot_big_pic                     = 1; %plot BIGPIC
computation_ROC          = 1; % compute ROC stat

source_clusters_err             = '';
draw_and_save_plots_err   = '';
draw_and_save_plots2_err = '';
plot_all_manual_spikes_err = '';
plot_big_pic_err                   = '';
computation_ROC_err        = '';
manual_spikes_nearest_spc_ica_err = '';
main190904_reduction_in_visualcluster_err = '';

%         % subj info
%         cortex           = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat, '\tess_cortex_pial_low.mat']));
%         MRI              = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat, '\subjectimage_T1.mat']));
%         Data             = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\', file_name,'\data_block001.mat']));
%         Data2           = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\B1C2_ii_run1_raw\data_block001.mat']));
%         channels      = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\channel_vectorview306_acc1.mat']));
%         G3                = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\headmodel_surf_os_meg.mat']));
%
% subj info

protocol_dir = 'E:\Valerii\Alex tes data 191104\EpilepsySolo'
cases_unique_for_anat = 'Sysoeva'

cortex           = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\anat\Sysoeva\tess_cortex_pial_low.mat');
MRI              = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\anat\Sysoeva\subjectimage_T1.mat');
Data             = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\data\Sysoeva\sleep-2_tsss\data_block001.mat');
%         Data2           = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\data\Sysoeva\sleep-2\data_block001.mat');
channels      = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\data\Sysoeva\@rawsleep-2\channel_vectorview306_acc1.mat');
G3                = load( 'E:\Valerii\Alex tes data 191104\EpilepsySolo\data\Sysoeva\@rawsleep-2_tsss\headmodel_surf_os_meg.mat');
labels           = extractfield(channels.Channel,'Name'); labels = labels(1:306)';
offset_time = Data.Time(1); %ms
 load('E:\Valerii\Alex tes data 191104\events_sleep-2_tsss.mat')

resultsdir_root = 'E:\Valerii\Alex tes data 191104\EpilepsySolo\Results\'
mkdir(resultsdir_root)
mkdir([resultsdir_root 'PICS'])
cases_unique_for_anat = 'Sysoeva'

Data.Events = events;
newdataset  = 0

 


 
 

%%
channel_type_loop = 2
switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
    case 1, channel_type = 'mag';  
                channel_idx     = 3:3:306;
    case 2, channel_type = 'grad'; 
                channel_idx     = setdiff(1:306, 3:3:306);
end
labels = labels(channel_idx) ;
 
%% 2. Spike detection

spikes_detection = 2,%detection_type%1:3

spikes_extraction = 'ICA_based';
decision = 0.9; % the amplitude threshold for decision (spike detection in ICA components)
f_low      = 3;    % highpass filter before the ICA decomposition
f_high     = 70; % lowpass filter before the ICA decomposition

ICA_spikes_mat = strcat([resultsdir_root, ...
    cases_unique_for_anat, '_ICA_based_', channel_type,'_EpiSolo.mat']);

if newdataset %which in this case == 1
    [spike_ind, picked_components, picked_comp_top, component_indicatior] = ...
        ICA_detection(Data, G3, channel_type, decision, f_low, f_high);
    load train;
    p = audioplayer(y, Fs);
    play(p, [1 (get(p, 'SampleRate') * 3)]);
    save(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
end

load(ICA_spikes_mat,'spike_ind', 'picked_components', 'picked_comp_top','component_indicatior')
% spike_ind = spike_ind(spike_inrd<600000-30 & spike_ind>41);
spike_clust = zeros(size(spike_ind))';
close all

ica_topo_time_detections_file = 
 ICA_plot_TOPO_Time_detections
 saveas(ica_topo_time_detections_file)

% 3. RAP-MUSIC (2) dipole fitting
if computation_source | newdataset
    
    if 1
        f_low_RAP  = 10;
        f_high_RAP = 200;
        spikydata = 0;
        %RAP = 'RAP'; corr_thresh = 0.99;
        RAP = 'not';
        if spikes_detection == 1
            corr_thresh = 0.0;
        else
            corr_thresh = CORR_THR; %0.95
        end
        % corr_thresh = back to quantile(ValMax, 0.95)
        [IndMax, ValMax, ind_m, spikeind] = spike_localization(spike_ind, Data, G3, ...
            channel_type, f_low_RAP, f_high_RAP, spikydata, picked_components, ...
            picked_comp_top, corr_thresh, RAP);
        xlim([0.95 1])
        
        spikes_fitted = [spikeind' spikeind'  ValMax(ind_m)' IndMax(ind_m)']%  [record time in s for brainstorm ¦ idx detecton sample  ¦ subspace corr ¦ leadfield_orig_index]
        spikes_fitted = sort(spikes_fitted,1)
        spikes_fitted(:,1) = spikes_fitted(:,1)/1000+offset_time
        
        plot_spikes_ER_TOPO(spikes_fitted, ([resultsdir_root 'PICS\']), Data, channel_type, f_low, f_high,  G3,  channels)
        
        [resultsdir_root,  '\sources_'  spikes_extraction '_' channel_type '.mat']
        %                         E:\Valerii\45_cases\B1C2\ASPIRE\CORR_THR_0.955\\sources_ICA_based_grad.mat
        
        save([resultsdir_root,  '\sources_'  spikes_extraction '_' channel_type '.mat'],...
            'IndMax','ValMax','ind_m','spikeind','spikes_fitted')
        
    else
        
        load([resultsdir_root,  '\sources_'  spikes_extraction '_' channel_type '.mat'],...
            'IndMax','ValMax','ind_m','spikeind')
        
    end
    
    %% clustering5
    Nmin = N_MIN;
    thr_dist = THR_DIST;
    
    Nmin = 5;
    thr_dist = 0.02;
    
    cluster = clustering(spike_ind, G3, Nmin, ValMax, IndMax, ind_m, ...
                            thr_dist, 1,cortex, RAP, spikeind, spike_clust);
    
    %% source space estimation
    f_low  = 3;
    f_high = 50;
    [spike_trials, maxamp_spike_av, channels_maxamp, spike_ts] = ...
        source_reconstruction(Data, G3, channel_type, cluster, ...
        f_low, f_high);
                    
     %% cluster plots
     f_low_vis  = 2; % bandpass filter for visualization
     f_high_vis = 50;
       save_param.resultsdir_root    =  resultsdir_root;
                    save_param.subj_name          = '';
                    save_param.results_subfolder  = [results_subfolder '\clusters'];
                    mkdir([ save_param.resultsdir_root save_param.results_subfolder])
                    save_param.spikes_extraction  = spikes_extraction;
     plot_clusters_191015(Data, channel_type, f_low_vis, f_high_vis, cortex, ...
                        spike_trials, maxamp_spike_av, spike_ts, cluster, ...
                        channels, G3, MRI, corr_thresh, save_param ) 
                    
                    
     %% 8. write clusters in csv to file
            if computation_clusters
                cluster_out  = cluster_out_spycirc(cluster, G3);
%                 csvwrite([resultsdir_root, subj_name, results_subfolder '\cluster_out_' spikes_extraction '_' channel_type '.csv'],cluster_out);
                csvwrite([resultsdir_root, '', results_subfolder '\cluster_out_' spikes_extraction '_' channel_type '.csv'],cluster_out);
                
                events = cluster_out(:,1)/1000;
                save(strcat([resultsdir_root, subj_name, '\Spyking_circus\', ...'_'
                    cases_unique_for_anat, '_', file_name_short, '\Aspire\',...
                    'EVENTS3_' spikes_extraction '_' file_name_short, '_', channel_type, '.mat']),'events');
                
            end           
            
            
            %% plot spikes of same cluster
            for cl = 1:length(cluster)
                
                spikes_fitted = [cluster_out(:,1) cluster_out(:,1) cluster_out(:,3) cluster_out(:,7)  cluster_out(:,2)]%  [record time in s for brainstorm ¦ idx detecton sample  ¦ subspace corr ¦ leadfield_orig_index ¦ cluster]
                spikes_fitted = sort(spikes_fitted,1)
                spikes_fitted(:,1) = spikes_fitted(:,1)/1000+offset_time
                
                spikes_fitted = spikes_fitted(find(spikes_fitted(:,end) == cl),:);
                
                mkdir([resultsdir_root 'PICS\cluster' num2str(cl) '\'])
                TOPO = plot_spikes_ER_TOPO(spikes_fitted, ([resultsdir_root 'PICS\cluster' num2str(cl) '\']), Data, channel_type, f_low, f_high,  G3,  channels)
                
                events = spikes_fitted(:,1)
                save ([resultsdir_root 'EVENTScluster' num2str(cl) ],'events')
                
            end
            
            
            
    %%
    figure(100)
    
    switch channel_type_loop % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'mag';
            subplot(221),
                histogram(ValMax),title('MAG'),ylabel('t__SSS'),
                xlim([0.8 1])
            subplot(223),
                histogram(ValMax2),ylabel('RAW'),
                xlim([0.8 1])
            
        case 2, channel_type = 'grad';
            subplot(222),
                histogram(ValMax),title('GRAD'),
                xlim([0.8 1])
            subplot(224),
                histogram(ValMax2),
                xlim([0.8 1])
            
    end
end



