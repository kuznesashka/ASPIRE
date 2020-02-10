function main()

%% -------------------------------------------------------------------------
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
paths.brainstorm = ;
paths.fieldtrip  = ;
addpath(genpath(paths.brainstorm));
addpath((paths.fieldtrip))
addpath(([paths.fieldtrip filesep 'plotting']));
addpath(([paths.fieldtrip filesep 'utilities' filesep 'private']));
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% -------------------------------------------------------------------------
%% Paths -- subject info
% -------------------------------------------------------------------------
% !!! File separator for the current platform
paths.anat 		 = ;
paths.data 		 = ;
paths.root 		 = ;

%settings_tommaso
paths.data_block =  ;
paths.subj_name  =  ;
paths.case 		 =  ;
paths.fname 	 =  ;
paths.sh_fname   =  ;


% subj info
paths.cortex 	 =  ;
paths.MRI 		 =  ;
paths.Data 		 =  ;
paths.channels   =  ;
paths.G3 		 =  ;

% paths to detections
paths.detections = ;

% -------------------------------------------------------------------------
%% Paths for saving


paths.sources_saving_path  =  ;
paths.clusters_saving_path =  ;
paths.results_saving_path  =  ;
paths.voxels_saving_path   =  ;
paths.affine_saving_path   =  ;

% -------------------------------------------------------------------------
%% Parameters main
% -------------------------------------------------------------------------
parameters.computation_source    = 1; % compute dipoles
parameters.computation_clusters  = 1; % compute and save clusters
parameters.run_ica               = 0; % not plot pictures


parameters.detection_type = [3]; %1-visual, 2-ICA, 3-SPC
parameters.channel_types = [2]; %1 - 'mag', 2 - 'grad'

% -------------------------------------------------------------------------
%% Parameters for detection
% -------------------------------------------------------------------------
parameters.detection.ICA.decision = 0.9; % the amplitude threshold for decision
parameters.detection.ICA.f_low 	  = 3; % bandpass filter before the ICA decomposition
parameters.detection.ICA.f_high   = 70;

% -------------------------------------------------------------------------
%% Parameters for RAP-MUSIC
% -------------------------------------------------------------------------
parameters.rap_music.f_low_RAP  = 10;
parameters.rap_music.f_high_RAP = 200;
parameters.rap_music.spikydata  = 0; % spikydata -- indicatior, showing whether you want to fit
parameters.rap_music.RAP        = 'not';
parameters.prctile              = 85; % prctile(ValMax,85); -- threshold for ICA and Spyking Circus
parameters.corr_thresh          = 0.95; % threshold for visual detections

% -------------------------------------------------------------------------
%% Parameters for clustering
% -------------------------------------------------------------------------
% THR_DIST - maximal distance from the center of the cluster (radius) in m
parameters.clustering.THR_DIST = 0.01;
% N_MIN - minimum number of sources in one cluster
parameters.clustering.N_MIN    = 3;

cortex          = load(paths.cortex);
MRI             = load(paths.MRI);
Data            = load(paths.Data);
channels        = load(paths.channels);
G3              = load(paths.G3);

%%
Voxels = cs_convert(MRI, 'scs', 'voxel', cortex.Vertices);
affine = MRI.InitTransf{2};

save(paths.voxels_saving_path, 'Voxels')
save(paths.affine_saving_path, 'affine')
%% run the main function
main_one_subject_run_from_python(cortex, Data, G3, MRI, channels, paths, parameters, channel_type, spikes_detection);

end

