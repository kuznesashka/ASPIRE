% 1. Export everything from brainstorm
subj_name = 'B1C2';
protocol_dir = '/home/ksasha/Documents/brainstorm_db/MEG_Tommaso/';
file_name = 'B1C2_ii_run1_raw_tsss_mc_art_corr';
channel_type = 'mag'; % channels you want to analyse ('grad' or 'mag') 

cortex = load(strcat([protocol_dir, 'anat/', subj_name,'/tess_cortex_pial_low.mat']));
MRI = load(strcat([protocol_dir, 'anat/', subj_name, '/subjectimage_T1.mat']));
Data = load(strcat([protocol_dir, 'data/', subj_name, '/', file_name,'/data_block001.mat']));
channels = load(strcat([protocol_dir, 'data/', subj_name, '/', file_name, '/channel_vectorview306_acc1.mat']));
G3 = load(strcat([protocol_dir, 'data/', subj_name, '/', file_name, '/headmodel_surf_os_meg.mat']));

% Here if you want to fit dipoles to the time-stamps computed with spiking
% circus, you just skip the step (2) and go to (3)

% 2. ICA-based spike detection
decision = 0.9; % the amplitude threshold for decision 
f_low = 3; % bandpass filter before the ICA decomposition
f_high = 70;
[spike_ind, picked_components, picked_comp_top] = ...
    ICA_detection(Data, G3, channel_type, decision, f_low, f_high);

% Upload the timestamps from spiking circus
% spcirc_data = csvread('Templates_B1C2.csv');
% spike_ind = spcirc_data(:,1);
% spike_ind = spike_ind(spike_ind<600000);

% 3. RAP-MUSIC dipole fitting
f_low = 10;
f_high = 200;
spikydata = 0;
[IndMax, ValMax, corr_thresh, ind_m] = spike_localization(spike_ind, Data, G3, channel_type, ...
    f_low, f_high, spikydata);

% 4. Clustering
thr_dist = 0.01; % maximal distance from the center of the cluster (radius)
Nmin = 8; % minimum number of sources in one cluster
cluster = clustering(spike_ind, G3, Nmin, ValMax, IndMax, ind_m, thr_dist, 1, cortex);

% 5. Activation on sources (broken part, but still you need to run it before the next plot)
f_low = 3;
f_high = 50;
[bf_ts, spike_sources, spike_sources_av] = beamformer_for_clusters(Data, ...
    G3, cluster, channel_type, 0);

% 6. Big plot
f_low = 2; % bandpass filter for visualization
f_high = 50;
epi_plot2(Data, channel_type, f_low, f_high, cortex, cluster, spike_ind, ...
    picked_components, picked_comp_top, spike_sources, bf_ts, MRI, G3, channels, corr_thresh, 1)
