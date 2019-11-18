% -------------------------------------------------------------------------
% Main with parameters
% -------------------------------------------------------------------------
% INPUTS:
%
% OUTPUTS:
%
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com



%load cases_files_190921.mat
%when two or more recordings 'B4Z2\Rec_01'
%paths.subj_name = cases_files.cases{case_n};
%paths.cases_unique_for_anat = cases_files.cases_unique{case_n};
%artifacts corrected file name
%paths.file_name = cases_files.file_names{case_n};
%short file name
%paths.file_name_short = cases_files.file_names_short{case_n};

%Debugging 
%dbstop if warning

%% Add libraries  
addpath(genpath('/Users/valery/MEG/brainstorm3'));
addpath /Users/valery/MEG/fieldtrip-20191008
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% -------------------------------------------------------------------------
%% Paths -- subject info
% -------------------------------------------------------------------------
% !!! File separator for the current platform
paths.anat = ['/Users/valery/MEG/EPILEPSY' filesep, ...
                'MEG_Tommaso' filesep 'anat' filesep]; % anatomy from Brainstorm
paths.data = ['/Users/valery/MEG/EPILEPSY' filesep, ...
                'MEG_Tommaso' filesep 'data' filesep]; % data from Brainstorm
paths.root = '/Users/valery/MEG/Cases/'; %45 cases folder

paths.subj_name = 'B1C2'; %when two or more recordings 'B4Z2\Rec_01'
paths.case = 'B1C2';
paths.fname = 'B1C2_ii_run1_raw_tsss_mc_art_corr'; %artifact corrected file name
paths.sh_fname = 'B1C2_ii_run1_raw_tsss_mc'; %short file name

mkdir([paths.root paths.subj_name filesep 'ASPIRE'])
mkdir([paths.root paths.subj_name filesep 'ASPIRE' filesep 'plots'])
mkdir([paths.root paths.subj_name filesep 'ASPIRE' filesep 'results'])

paths.detections = [paths.root paths.subj_name filesep 'ASPIRE', ...
                    filesep 'detections' filesep];
paths.plots = [paths.root paths.subj_name filesep 'ASPIRE' filesep, ...
                    'plots' filesep];
paths.results = [paths.root paths.subj_name filesep 'ASPIRE' filesep, ...
                    'results' filesep];

% subj info
paths.cortex = strcat([paths.anat paths.case filesep 'tess_cortex_pial_low.mat']);
paths.MRI = strcat([paths.anat paths.case filesep 'subjectimage_T1.mat']);
paths.Data = strcat([paths.data paths.case filesep paths.fname, ...
                    filesep 'data_block001.mat']);
paths.channels = strcat([paths.data paths.case filesep '@default_study', ...
                         filesep 'channel_vectorview306_acc1.mat']);
paths.G3 = strcat([paths.data paths.case filesep '@default_study', ...
                   filesep 'headmodel_surf_os_meg.mat']);

% paths to detections
% path_vis_detections -- path to the visual detections file (csv)
paths.path_vis_detections = [paths.detections 'Manual_spikes_' paths.sh_fname '.csv'];
% path_ICA_detections -- path to the ICA detections (mat)
paths.path_ICA_detections = [paths.detections 'ICA_detections_' paths.sh_fname];
% path_SPC_detections -- path to the Spyking Circus detections (csv)
paths.path_SPC_detections = [paths.detections 'Templates_' paths.sh_fname];

% -------------------------------------------------------------------------
%% Paths for saving
% -------------------------------------------------------------------------
% Path for sources saving without [spikes_extraction '_' channel_type '.mat']
paths.sources_saving_path = [paths.results 'sources_'];
% Path for clusters saving without [spikes_extraction '_' channel_type '.csv']
paths.path_cluster_out = [paths.results 'cluster_out_'];
% Path for results saving without [spikes_extraction '_' channel_type '.mat']
paths.results_saving_path = [paths.results 'results_'];

% Save clusters plots
mkdir([paths.root paths.subj_name filesep 'ASPIRE' filesep 'plots' filesep 'clusters'])
paths.save_cluster_plots = [paths.plots filesep 'clusters' filesep ];

% ROC saving path
mkdir([paths.root 'ROC'])
paths.roc = [paths.root 'ROC' filesep];
paths.roc_xlsx_fname = [paths.roc 'ROC.xlsx'];
paths.roc_labels_xlsx_fname  = [paths.roc 'Labels.xlsx'];

% Path for saving big picture
paths.bigpic_saving_path = [paths.plots paths.case '.bmp'];
paths.overlap_saving_path = [paths.plots 'overlap_'];

% -------------------------------------------------------------------------
%% Parameters main
% -------------------------------------------------------------------------
parameters.computation_source    = 1; % compute dipoles
parameters.computation_clusters  = 1; % compute and save clusters
parameters.draw_and_save_plots   = 1; % plot clusters
parameters.save_results			 = 1; % save results in file
parameters.computation_ROC       = 1; % compute ROC stat
parameters.plot_big_pic          = 0; % Plot all detections on the big plot
parameters.mute_mode             = 1; % if 0 - plot clickable clusters plot 
% in the not mute mode program will pause until a figure close
parameters.newdataset            = 0; % not plot pictures
parameters.propagation_probability = 0; % all to all cluster propagation probability
parameters.compute_overlap       = 1; % Overlap between detections


parameters.detection_type = [1 2 3]; %1-visual, 2-ICA, 3-SPC

% -------------------------------------------------------------------------
%% Parameters for detection
% -------------------------------------------------------------------------
parameters.detection.ICA.spikes_extraction = 'ICA_based';
parameters.detection.ICA.decision = 0.9; % the amplitude threshold for decision
parameters.detection.ICA.f_low = 3; % bandpass filter before the ICA decomposition
parameters.detection.ICA.f_high = 70;

parameters.detection.SPC.spikes_extraction = 'SpyCir_based';
parameters.detection.visual.spikes_extraction = 'visual';

% -------------------------------------------------------------------------
%% Parameters for RAP-MUSIC
% -------------------------------------------------------------------------
parameters.rap_music.f_low_RAP  = 10;
parameters.rap_music.f_high_RAP = 200;
parameters.rap_music.spikydata = 0; % spikydata -- indicatior, showing whether you want to fit
parameters.rap_music.RAP = 'not';
parameters.prctile = 85; % prctile(ValMax,85); -- threshold for ICA and Spyking Circus
parameters.corr_thresh = 0.95; % threshold for visual detections

% -------------------------------------------------------------------------
%% Parameters for clustering
% -------------------------------------------------------------------------
% THR_DIST - maximal distance from the center of the cluster (radius) in m
parameters.clustering.THR_DIST =  0.01;
% N_MIN - minimum number of sources in one cluster
parameters.clustering.N_MIN = 3;

%% Parameters for plots
parameters.draw.f_low  = 3;
parameters.draw.f_high = 50;
parameters.draw.f_low_vis  = 2;
parameters.draw.f_high_vis = 50;
parameters.draw.save_clusters = 0; %save all cluster as separate files


cortex          = load(paths.cortex);
MRI             = load(paths.MRI);
Data            = load(paths.Data);
channels        = load(paths.channels);
G3              = load(paths.G3);

main_one_subject(cortex, Data, G3, MRI, channels, paths, parameters)




