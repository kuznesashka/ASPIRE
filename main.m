% -------------------------------------------------------------------------
% Main with parameters
% -------------------------------------------------------------------------
% INPUTS:
%   cases_files_190921.mat -- all cases with path, subject name etc.
%   param_matrix190809.mat -- parametrs for iteration
%   newdataset            = 0;
%   computation_source    = 1; % compute dipoles
%   computation_clusters  = 1; % compute clustering
%   draw_and_save_plots   = 0; % plot clusters
%   draw_and_save_plots2  = 0; % save 
%   computation_ROC       = 1; % compute ROC stat
%
%
% OUTPUTS:
%
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com



load cases_files_190921.mat


%% Paths -- subject info
paths.hdisk = 'D:\';
paths.subj_name = cases_files.cases{case_n};
paths.cases_unique_for_anat = cases_files.cases_unique{case_n};
paths.protocol_dir = [hdisk 'Valerii\EPILEPSY\MEG_Tommaso\'];
paths.file_name = cases_files.file_names{case_n};
paths.file_name_short = cases_files.file_names_short{case_n};
paths.resultsdir_root = [hdisk 'Valerii\45_cases\'];
paths.results_subfolder = '\ASPIRE\';
mkdir([paths.resultsdir_root paths.subj_name '\' paths.results_subfolder])


% paths to detections
% path_vis_detections -- path to the visual detections file (csv)
paths.path_vis_detections = 
% path_ICA_detections -- path to the ICA detections (mat)
paths.path_ICA_detections = 
% path_SPC_detections -- path to the Spyking Circus detections (csv)
paths.path_SPC_detections = 

% subj info
paths.cortex = strcat([protocol_dir, 'anat\', cases_unique_for_anat, ...
                            '\tess_cortex_pial_low.mat'])
paths.MRI = strcat([protocol_dir, 'anat\', cases_unique_for_anat, ...
                            '\subjectimage_T1.mat'])
paths.Data = strcat([protocol_dir, 'data\', cases_unique_for_anat, ...
                            '\', file_name,'\data_block001.mat'])
paths.channels = strcat([protocol_dir, 'data\', cases_unique_for_anat, ...
                            '\@default_study', '\channel_vectorview306_acc1.mat'])
paths.G3 = strcat([protocol_dir, 'data\', cases_unique_for_anat, ...
                            '\@default_study', '\headmodel_surf_os_meg.mat'])


%% Paths for saving
% Path for sources saving without [spikes_extraction '_' channel_type '.mat']
paths.sources_saving_path = [resultsdir_root subj_name results_subfolder '\sources_'];
% Path for clusters saving without [spikes_extraction '_' channel_type '.csv']
paths.path_cluster_out = [resultsdir_root subj_name results_subfolder '\cluster_out_'];
% Path for results saving without [spikes_extraction '_' channel_type '.mat']
paths.results_saving_path = [resultsdir_root subj_name results_subfolder '\results_'];

% ROC saving path
paths.roc_xlsx_fname = [resultsdir_root 'Aspire_ROC\ROC_COR_TR_' num2str(CORR_THR) '.xlsx'];
paths.roc_labels_xlsx_fname  = [resultsdir_root 'Aspire_ROC\Labels_COR_TR_' num2str(CORR_THR) '.xlsx'];
% Path for saving big picture
paths.bigpic_saving_path = [resultsdir_root, subj_name, '\ASPIRE\CORR_THR_', ...
    num2str(CORR_THR), '_', xlsTAB '.bmp'] 


% Path for sources saving without [spikes_extraction '_' channel_type '.mat']
paths.sources_saving_path = [resultsdir_root subj_name results_subfolder '\sources_'];
% Path for clusters saving without [spikes_extraction '_' channel_type '.csv']
paths.path_cluster_out = [resultsdir_root subj_name results_subfolder '\cluster_out_'];
% Path for results saving without [spikes_extraction '_' channel_type '.mat']
paths.results_saving_path = [resultsdir_root subj_name results_subfolder '\results_'];


%% Parameters main
parameters.computation_source    = 1; % compute dipoles
parameters.computation_clusters  = 1; % compute clustering
parameters.draw_and_save_plots   = 0; % plot clusters
parameters.draw_and_save_plots2  = 0; % save clustering
parameters.computation_ROC       = 1; % compute ROC stat
parameters.plot_big_pic      = 1; % 
parameters.mute_mode      = 1; % not plot pictures
parameters.newdataset     = 1; % not plot pictures

cortex          = load(paths.cortex);
MRI             = load(paths.MRI);
Data            = load(paths.Data);
channels        = load(paths.channels);
G3              = load(paths.G3);

parameters.detection_type = [1 2 3];

%% Parameters detection
parameters.detection.ICA.spikes_extraction = 'ICA_based';
parameters.detection.ICA.decision = 0.9; % the amplitude threshold for decision
parameters.detection.ICA.f_low = 3; % bandpass filter before the ICA decomposition
parameters.detection.ICA.f_high = 70;

parameters.detection.SPC.spikes_extraction = 'SpyCir_based';
parameters.detection.visual.spikes_extraction = 'visual';

%% Parameters RAP-MUSIC
parameters.rap_music.f_low_RAP  = 10;
parameters.rap_music.f_high_RAP = 200;
parameters.rap_music.spikydata = 0; % spikydata -- indicatior, showing whether you want to fit
parameters.rap_music.RAP = 'not';


%% Parameters clustering
% THR_DIST - maximal distance from the center of the cluster (radius) in m
parameters.clustering.THR_DIST =  0.01;
% N_MIN - minimum number of sources in one cluster
parameters.N_MIN = 3;
parameters.CORR_THR = 0.95;


%% Parameters plots
parameters.draw.f_low  = 3;
parameters.draw.f_high = 50;
parameters.draw.f_low_vis  = 2;
parameters.draw.f_high_vis = 50;






main_one_subject(cortex, Data, G3, paths, parameters)
















