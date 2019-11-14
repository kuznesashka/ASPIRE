function ROC(subj_name, detection_type, path_cluster_out, ...
             cortex, roc_xlsx_fname, roc_labels_xlsx_fname)
%
% -------------------------------------------------------------------------
% Computation and saving ROC
% -------------------------------------------------------------------------
% ICA_grad
% ICA_mag
% SPC_grad
% SPC_mag
%
% roc_xlsx_fname - output file for ROC tables
% roc_labels_xlsx_fname - output file for labels tables
%
% OUTPUTS:
%
%_______________________________________________________
%
%
% Load Spyking Circus detected timestamps
SPC_grad = load([path_cluster_out 'SpyCir_based_grad.csv']);
SPC_mag = load([path_cluster_out 'SpyCir_based_mag.csv']);

% Load ICA detected timestamps
ICA_grad = load([path_cluster_out 'ICA_based_grad.csv']);
ICA_mag = load([path_cluster_out 'ICA_based_mag.csv']);

% Load visual timestamps
visual = load([path_cluster_out 'visual_grad.csv']);
    
xlsTAB = subj_name;
backsalsh = find(xlsTAB=='\');
if backsalsh
    xlsTAB(backsalsh)='_';
end


for spikes_detection = detection_type%1:3
    
    switch spikes_detection
        
        case 2
            
            % ICA
            fname = roc_xlsx_fname;
            fname_labels = roc_labels_xlsx_fname;
            
            ICA_visual = visual;
            ICA_visual_mag = visual;
            ICA_visual_grad = visual;
            
            % ROC in time domain
            ICA_roc = roc_time(ICA_visual(ICA_visual(:,1)<600000,:), ICA_grad, ICA_mag);
            ICA_roc_T = cell2table(ICA_roc(2:end,:));
            ICA_roc_T.Properties.VariableNames = ICA_roc(1,:);
            ICA_roc_T.Properties.VariableNames(1) = {'ICA'};

            % ROC in spatial domain
            ICA_roc_spatial = roc_spatial(ICA_visual_grad, ICA_visual_mag, ICA_grad, ICA_mag);
            ICA_roc_spatial_T = cell2table(ICA_roc_spatial(2:end,:));
            ICA_roc_spatial_T.Properties.VariableNames = ICA_roc_spatial(1,:);
            ICA_roc_spatial_T.Properties.VariableNames(1) = {'ICA'};

            % ROC in anatomical labels domain
            [ICA_labels_results, ICA_roc_labels] = roc_labels(ICA_visual_mag, ICA_visual_grad, ICA_mag, ...
                ICA_grad, cortex.Atlas);

            ICA_labels_results_T = cell2table(ICA_labels_results(2:end,:));
            ICA_labels_results_T.Properties.VariableNames = ICA_labels_results(1,:);
            ICA_labels_results_T.Properties.VariableNames(1) = {'ICA'};

            ICA_roc_labels_T = cell2table(ICA_roc_labels(2:end,:));
            ICA_roc_labels_T.Properties.VariableNames = ICA_roc_labels(1,:);
            ICA_roc_labels_T.Properties.VariableNames(1) = {'ICA'};

            % Save results in the Excel file
            writetable(ICA_roc_T,fname,'Sheet',xlsTAB,'Range','A6');
            writetable(ICA_roc_spatial_T,fname,'Sheet',xlsTAB,'Range','A16');
            writetable(ICA_roc_labels_T,fname,'Sheet',xlsTAB,'Range','A26');
            writetable(ICA_labels_results_T,fname_labels,'Sheet',xlsTAB,'Range','I1');

            % Only for Windows
            %xlswrite(fname,ICA_roc,xlsTAB,'A6');
            %xlswrite(fname,{'ICA'},xlsTAB,'A6');
            %xlswrite(fname,ICA_roc_spatial,xlsTAB,'A16');
            %xlswrite(fname,{'ICA'},xlsTAB,'A16');
            %xlswrite(fname,ICA_roc_labels,xlsTAB,'A26');
            %xlswrite(fname,{'ICA'},xlsTAB,'A26');
            %xlswrite(fname_labels,ICA_labels_results,xlsTAB,'I1');
            %xlswrite(fname_labels,{'ICA'},xlsTAB,'I1');
            
        case 3
            % Spyking Circus
            fname = roc_xlsx_fname;
            fname_labels = roc_labels_xlsx_fname;
            
            SPC_visual = visual;
            SPC_visual_mag = visual;
            SPC_visual_grad = visual;
            
            % ROC in time domain
            SPC_roc = roc_time(SPC_visual(SPC_visual(:,1)<600000,:), SPC_grad, SPC_mag);
            SPC_roc_T = cell2table(SPC_roc(2:end,:));
            SPC_roc_T.Properties.VariableNames = SPC_roc(1,:);
            SPC_roc_T.Properties.VariableNames(1) = {'SPC'};

            % ROC in spatial domain
            SPC_roc_spatial = roc_spatial(SPC_visual_grad, SPC_visual_mag, SPC_grad, SPC_mag);
            SPC_roc_spatial_T = cell2table(SPC_roc_spatial(2:end,:));
            SPC_roc_spatial_T.Properties.VariableNames = SPC_roc_spatial(1,:);
            SPC_roc_spatial_T.Properties.VariableNames(1) = {'SPC'};

            % ROC in anatomical labels domain
            [SPC_labels_results, SPC_roc_labels] = roc_labels(SPC_visual_mag, SPC_visual_grad, SPC_mag, ...
                SPC_grad, cortex.Atlas);

            SPC_labels_results_T = cell2table(SPC_labels_results(2:end,:));
            SPC_labels_results_T.Properties.VariableNames = SPC_labels_results(1,:);
            SPC_labels_results_T.Properties.VariableNames(1) = {'SPC'};

            SPC_roc_labels_T = cell2table(SPC_roc_labels(2:end,:));
            SPC_roc_labels_T.Properties.VariableNames = SPC_roc_labels(1,:);
            SPC_roc_labels_T.Properties.VariableNames(1) = {'SPC'};

            % Save results in the Excel file
            writetable(SPC_roc_T,fname,'Sheet',xlsTAB,'Range','A1');
            writetable(SPC_roc_spatial_T,fname,'Sheet',xlsTAB,'Range','A11');
            writetable(SPC_roc_labels_T,fname,'Sheet',xlsTAB,'Range','A21');
            writetable(SPC_labels_results_T,fname_labels,'Sheet',xlsTAB,'Range','A1');
            
            % Only for Windows
            %xlswrite(fname,SPC_roc,xlsTAB,'A1');
            %xlswrite(fname,{'SPC'},xlsTAB,'A1');
            %xlswrite(fname,SPC_roc_spatial,xlsTAB,'A11');
            %xlswrite(fname,{'SPC'},xlsTAB,'A11');
            %xlswrite(fname,SPC_roc_labels,xlsTAB,'A21');
            %xlswrite(fname,{'SPC'},xlsTAB,'A21');
            %xlswrite(fname_labels,SPC_labels_results,xlsTAB,'A1');
            %xlswrite(fname_labels,{'SPC'},xlsTAB,'A1');
            
    end
end

end