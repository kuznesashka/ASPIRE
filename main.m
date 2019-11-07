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
%   draw_and_save_plots2  = 0; % save clustering
%   computation_ROC       = 1; % compute ROC stat
%
%
% OUTPUTS:
%
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com



load cases_files_190921.mat

hdisk = 'D:\';
subj_name = cases_files.cases{case_n};
cases_unique_for_anat = cases_files.cases_unique{case_n};
protocol_dir = [hdisk 'Valerii\EPILEPSY\MEG_Tommaso\'];
file_name = cases_files.file_names{case_n};
file_name_short = cases_files.file_names_short{case_n};
resultsdir_root = [hdisk 'Valerii\45_cases\'];
results_subfolder = '\ASPIRE\';
mkdir([resultsdir_root subj_name '\' results_subfolder])


computation_source    = 1; % compute dipoles
computation_clusters  = 1; % compute clustering
draw_and_save_plots   = 0; % plot clusters
draw_and_save_plots2  = 0; % save clustering
computation_ROC       = 1; % compute ROC stat

% subj info
cortex         = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat,'\tess_cortex_pial_low.mat']));
MRI             = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat, '\subjectimage_T1.mat']));
Data            = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\', file_name,'\data_block001.mat']));
channels        = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\channel_vectorview306_acc1.mat']));
G3              = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\headmodel_surf_os_meg.mat']));

detection_type = [1 2 3];


roc_xlsx_fname = [resultsdir_root 'Aspire_ROC\ROC_COR_TR_' num2str(CORR_THR) '.xlsx'];
roc_labels_xlsx_fname = [resultsdir_root 'Aspire_ROC\Labels_COR_TR_' num2str(CORR_THR) '.xlsx'];

main_one_subject()







