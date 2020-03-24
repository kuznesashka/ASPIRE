function main_run_from_python(paths_params)

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
load(paths_params)

%% Add libraries  
paths.brainstorm = paths_params.brainstorm;
paths.fieldtrip  = paths_params.fieldtrip;
addpath(genpath(paths.brainstorm));
addpath((paths.fieldtrip))
addpath(([paths.fieldtrip filesep 'plotting']));
addpath(([paths.fieldtrip filesep 'utilities' filesep 'private']));
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if paths_params.propagation > 0
    addpath(([paths.fieldtrip filesep 'utilities']));
end
% -------------------------------------------------------------------------
%% Paths -- subject info
% -------------------------------------------------------------------------

% subj info
paths.cortex 	 =  paths_params.cortex;
paths.MRI 		 =  paths_params.MRI;
paths.channels   =  paths_params.channels;
paths.G3 		 =  paths_params.G3;

% paths to detections
paths.detections = paths_params.detections;

% -------------------------------------------------------------------------
%% Paths for saving


paths.sources_saving_path  =  paths_params.sources_saving_path;
paths.clusters_saving_path =  paths_params.clusters_saving_path;
paths.parameters_saving_path  =  paths_params.parameters_saving_path;
paths.voxels_saving_path   =  paths_params.voxels_saving_path;
paths.affine_saving_path   =  paths_params.affine_saving_path;

% -------------------------------------------------------------------------
%% Parameters main
% -------------------------------------------------------------------------
parameters.N_data 				 = paths_params.N_data;
for data_n = 1:parameters.N_data
    switch data_n
        case 1, parameters.Data_0 = paths_params.Data_0;
        case 2, parameters.Data_1 = paths_params.Data_1;
        case 3, parameters.Data_2 = paths_params.Data_2;
    end
end

parameters.computation_source    = 1; % compute dipoles
parameters.computation_clusters  = 1; % compute and save clusters
parameters.run_ica               = 0; % not plot pictures


parameters.spikes_detection = [paths_params.detection_type]; %1-visual, 2-ICA, 3-SPC, 4-alphaCSC
parameters.channel_type  = [paths_params.channel_types]; %1 - 'mag', 2 - 'grad'

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

if paths_params.propagation == 0
    parameters.prctile              = 95; % prctile(ValMax,85); -- threshold for ICA and Spyking Circus
    parameters.corr_thresh          = 0.95; % threshold for visual detections
else
    parameters.prctile              = 85; % prctile(ValMax,85); -- threshold for ICA and Spyking Circus
    parameters.corr_thresh          = 0.8; % threshold for visual detections
end
% -------------------------------------------------------------------------
%% Parameters for clustering
% -------------------------------------------------------------------------
% THR_DIST - maximal distance from the center of the cluster (radius) in m
if paths_params.propagation == 0
	parameters.clustering.THR_DIST = 0.01;
else
	parameters.clustering.THR_DIST = 0.02;
end
% N_MIN - minimum number of sources in one cluster
parameters.clustering.N_MIN    = 5;

cortex          = load(paths.cortex);
MRI             = load(paths.MRI);
channels        = load(paths.channels);
G3              = load(paths.G3);

%%
Voxels = cs_convert(MRI, 'scs', 'voxel', cortex.Vertices);
Voxels_mni = cs_convert(MRI, 'scs', 'mni', cortex.Vertices);
affine_scs = [MRI.SCS.R MRI.SCS.T; 0 0 0 1];
affine_ncs = [MRI.NCS.R MRI.NCS.T; 0 0 0 1];

save(paths.voxels_saving_path, 'Voxels')
save(paths_params.voxels_mni_saving_path, 'Voxels_mni')
save(paths.affine_saving_path, 'affine_scs')
save(paths_params.affine_mni_saving_path, 'affine_ncs')
%% run the main function
if paths_params.propagation == 0
	main_one_subject_run_from_python(cortex, G3, MRI, channels, paths, parameters);
elseif paths_params.propagation == 1
	parameters.t1 = paths_params.t1;
	parameters.t2 = paths_params.t2;
	parameters.t3 = paths_params.t3;
    paths.sources_saving_path_t1 = paths_params.sources_saving_path_t1;
    paths.sources_saving_path_t3 = paths_params.sources_saving_path_t3;
	main_one_subject_propagation_run_from_python(cortex, G3, MRI, channels, paths, parameters);
elseif paths_params.propagation == 2 %dipole fit
    switch parameters.channel_type % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'grad';
            channel_idx     = setdiff(1:306, 3:3:306);
        case 2, channel_type = 'mag';
            channel_idx     = 3:3:306;
    end
    Data = load(parameters.Data_0);
    evoked_data = load(paths_params.evoked);
    t1= paths_params.t1;
    t2= paths_params.t2;
    t3= paths_params.t3;
    t4= paths_params.t4;
    t = [t1 int64((t1+t2)/2) t2 int64((t2+t3)/2) t3 int64((t3+t4)/2) t4];
    dip_fit_average(Data, evoked_data, G3, channels, channel_idx, t,...
                    MRI, cortex, paths_params.evoked_saving_path)
elseif paths_params.propagation == 3 %beamforming
    Data = load(paths_params.atoms_epoch_data);
    Data = Data.data;
    switch parameters.channel_type % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'grad';
            channel_idx     = setdiff(1:306, 3:3:306);
        case 2, channel_type = 'mag';
            channel_idx     = 3:3:306;
            % Data = Data * 100;
    end
    dip_ind = paths_params.dip_ind;
    VE   = source_reconstruction_atoms(Data', G3, channel_idx, dip_ind);
    save(paths_params.atoms_soure_data, 'VE')
%     chunk_len = length(VE)/99;
%     avg_VE = zeros(length(dip_ind), chunk_len);
%     for chunk = 1:length(VE)/chunk_len
%         chunk_begin = 1 + (chunk-1)*chunk_len;
%         chunk_end   = chunk_len + (chunk-1)*chunk_len;
%         avg_VE(:,:) = avg_VE(:,:) + VE(:,chunk_begin:chunk_end);
%     end
%     for i = 1:8
%         figure(i)
%         plot(avg_VE(i,:)/99);
%     end
%    VE   = source_reconstruction_atoms(data, G3, channel_idx, dip_ind, 1);

elseif paths_params.propagation == 4 %RAP MUSIC
    evoked_data = load(paths_params.evoked);
    switch parameters.channel_type % channels you want to analyse ('grad' or 'mag')
        case 1, channel_type = 'grad';
            channel_idx     = setdiff(1:306, 3:3:306);
        case 2, channel_type = 'mag';
            channel_idx     = 3:3:306;
            % Data = Data * 100;
    end
    [Valmax, Indmax, Sources] = RAP_MUSIC_scan_atoms(evoked_data.evoked, ...
        parameters.corr_thresh, G3, channel_idx);
    save(paths_params.rap_save, 'Valmax', 'Indmax', 'Sources')

end
end

