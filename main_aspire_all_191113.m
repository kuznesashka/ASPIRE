% -------------------------------------------------------------------------
% Main with parameters
% -------------------------------------------------------------------------
% INPUTS:
%   cases_files_190921.mat -- all cases with path, subject name etc.
%   param_matrix190809.mat -- parametrs for iteration
%   newdataset            = 0;
%   computation_source    = 1; % compute dipoles
%	computation_clusters  = 1; % compute clustering
%	draw_and_save_plots   = 0; % plot clusters
%	draw_and_save_plots2  = 0; % save clustering
%	computation_ROC       = 1; % compute ROC stat
%   plot_manual_spikes    = 1; % plot all manual spikes
%
%
% OUTPUTS:
%
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com




%% libraries
START_ASPIRE
%% 0. PAREMETERS

%load cases names
load cases_files_191022.mat

%!!IN CLUSTERIN draw = 0
% corr_thresh = back to quantile(ValMax, quant), but 300 spikes

% CORR_THR__ = [0.9 .925 .95]'
% THR_DIST__ = [0.01 0.02 0.03]'; % maximal distance from the center of the cluster (radius) in m
% N_MIN__    = [3 5 8]'; %
load param_matrix190809.mat

OUTPUT_visclust = []; % how many spikes survive after clustering
ROC_line_excel_spatial_ICA = [ ];
ROC_line_excel_labels_ICA  = [ ];
ROC_line_excel_spatial_SPC = [ ];
ROC_line_excel_labels_SPC  = [ ];
ROC_line_excel_time_ICA    = [ ];
ROC_line_excel_time_SPC    = [ ];

%% Parameters iteration
for case_n =  [2:4 7:length(cases_files.file_names_short)] %1:4
    
    newdataset = 1;
    %
    %         for case_n = [11:length(file_names_short)        ]
    
    STAT_ICA_temp = [];
    STAT_ICA_spac = [];
    STAT_ICA_labe = [];
    
    STAT_SPC_temp = [];
    STAT_SPC_spac = [];
    STAT_SPC_labe = [];
    
    for c_tr = [0.955] %var1 = 19 %:length(param_matrix)
        
        
        var1 = 19;
        CORR_THR = c_tr; %param_matrix(var1,1);
        THR_DIST = 0.02;%param_matrix(var1,2);
        N_MIN        = 5;%param_matrix(var1,3);
        
        
        close all
        %         try
        subj_name = cases_files.cases{case_n};
        cases_unique_for_anat = cases_files.cases_unique{case_n};
        %protocol_dir = '\\172.16.118.134\TFedele\Tommaso\Matlab Toolboxes\EPILEPSY\MEG_Tommaso\';
        protocol_dir = [hdisk 'Valerii\EPILEPSY\MEG_Tommaso\'];
        file_name = cases_files.file_names{case_n};
        file_name_short = cases_files.file_names_short{case_n};
        resultsdir_root = [hdisk 'Valerii\45_cases\'];
%         results_subfolder = ['\ASPIRE\CORR_THR_' num2str(CORR_THR) '\'];
%         mkdir([resultsdir_root subj_name results_subfolder])
        mkdir([resultsdir_root subj_name '\ASPIRE\'])
        
        computation_source    = 1; % compute dipoles
        computation_clusters  = 1; % compute clustering
        draw_and_save_plots   = 1; % plot clusters
        draw_and_save_plots2  = 1; % save clustering
        plot_all_manual_spikes    = 0; % plot all manual spikes
        plot_big_pic          =  0; %plot BIGPIC
        computation_ROC       = 0; % compute ROC stat
        
        source_clusters_err   = '';
        draw_and_save_plots_err = '';
        draw_and_save_plots2_err = '';
        plot_all_manual_spikes_err = '';
        plot_big_pic_err = '';
        computation_ROC_err = '';
        manual_spikes_nearest_spc_ica_err = '';
        main190904_reduction_in_visualcluster_err = '';
        
        
        % subj info
        cortex         = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat,'\tess_cortex_pial_low.mat']));
        MRI             = load(strcat([protocol_dir, 'anat\', cases_unique_for_anat, '\subjectimage_T1.mat']));
        Data           = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\', file_name,'\data_block001.mat']));
        channels        = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\channel_vectorview306_acc1.mat']));
        G3              = load(strcat([protocol_dir, 'data\', cases_unique_for_anat, '\@default_study', '\headmodel_surf_os_meg.mat']));
        labels           = extractfield(channels.Channel,'Name'); labels = labels(1:306)';
        offset_time = Data.Time(1); %ms
%         detection_type = [1 2 3];
                detection_type = [2];

        %% All steps for one case
%         main_one_run_191015
       main_one_run_191113
    
        
    
        
        %% Write errors
        xlsTAB = subj_name;
        backsalsh = find(xlsTAB=='\');
        if backsalsh
            xlsTAB(backsalsh)='_';
        end
        
        fname = [resultsdir_root 'Aspire_ROC\Errors_COR_TR_' num2str(CORR_THR) '.xlsx'];

        xlswrite(fname,{'manual_spikes_nearest_spc_ica'},xlsTAB,'A7');
        xlswrite(fname,{manual_spikes_nearest_spc_ica_err},xlsTAB,'D7');
        xlswrite(fname,{'main190904_reduction_in_visualcluster'},xlsTAB,'A8');
        xlswrite(fname,{main190904_reduction_in_visualcluster_err},xlsTAB,'D8');
        
    end
    
end



%% Write information about visual spikes
if 0
% write to E:\Valerii\45_cases\Aspire_ROC\ROCallsubjects190905.xlsx
% main190904_reduction_in_visualcluster
excelfileallsubjs = [resultsdir_root 'Aspire_ROC\ROCallsubjects190921.xlsx'];
file_names_short = cases_files.file_names_short;

file_name_short = cases_files.file_names_short{case_n};
xlswrite(excelfileallsubjs,ROC_line_excel_spatial_SPC,'Spatial','C19')
xlswrite(excelfileallsubjs,ROC_line_excel_spatial_ICA,'Spatial','C3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Spatial','B3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Spatial','B19')

xlswrite(excelfileallsubjs,ROC_line_excel_labels_SPC,'Labels','C19')
xlswrite(excelfileallsubjs,ROC_line_excel_labels_ICA,'Labels','C3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Labels','B3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Labels','B19')

xlswrite(excelfileallsubjs,ROC_line_excel_time_SPC,'Time','C19')
xlswrite(excelfileallsubjs,ROC_line_excel_time_ICA,'Time','C3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Time','B3')
xlswrite(excelfileallsubjs,file_names_short([1:4 7:13])','Time','B19')

xlswrite(excelfileallsubjs,OUTPUT_visclust,'SPIKES','B2')
xlswrite(excelfileallsubjs,{'subj','sens','vis','thr','vis>thr','N sp clustr','spikes in sp clustr','N label clustr','spikes in label clustr'},'SPIKES','B1')
xlswrite(excelfileallsubjs,reshape(repmat(file_names_short([1:4 7:13]),2,1),[],1),'SPIKES','A2')


figure('visible','off'),
subplot(211)
boxplot(ROC_line_excel_spatial_ICA(:,[1 3 5]))
subplot(212)
boxplot(ROC_line_excel_spatial_SPC(:,[1 3 5]))


reshape(repmat(file_names_short([1:4 7:13]),2,1),[],1)

end
