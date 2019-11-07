function ROC(detection_type, ICA_grad, ICA_mag, SPC_grad, SPC_mag, visual, ...
    cortex, roc_xlsx_fname, roc_labels_xlsx_fname, )
%
% -------------------------------------------------------------------------
% Computation and saving ROC
% -------------------------------------------------------------------------
% INPUTS:
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
            
            % ROC in spatial domain
            ICA_roc_spatial = roc_spatial(ICA_visual_grad, ICA_visual_mag, ICA_grad, ICA_mag);
            
            % ROC in anatomical labels domain
            [ICA_labels_results, ICA_roc_labels] = roc_labels(ICA_visual_mag, ICA_visual_grad, ICA_mag, ...
                ICA_grad, cortex.Atlas);
            
            % Save results in the Excel file
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
            fname = roc_xlsx_fname;
            fname_labels = roc_labels_xlsx_fname;
            
            SPC_visual = visual;
            SPC_visual_mag = visual;
            SPC_visual_grad = visual;
            
            % ROC in time domain
            SPC_roc = roc_time(SPC_visual(SPC_visual(:,1)<600000,:), SPC_grad, SPC_mag);
            
            % ROC in spatial domain
            SPC_roc_spatial = roc_spatial(SPC_visual_grad, SPC_visual_mag, SPC_grad, SPC_mag);
            
            % ROC in anatomical labels domain
            [SPC_labels_results, SPC_roc_labels] = roc_labels(SPC_visual_mag, SPC_visual_grad, SPC_mag, ...
                SPC_grad, cortex.Atlas);
            
            % Save results in the Excel file
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

end