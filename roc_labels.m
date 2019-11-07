function [labels_results, roc] = roc_curve_labels(v_mag, v_grad, mag, grad, atlas)

% -------------------------------------------------------------------------
% anatomical ROC curve
% -------------------------------------------------------------------------
% INPUTS:
%   v_mag, v_grad, mag, grad -- results with vertices (column 7)
%   atlas -- cortex.Atlas
%
% OUTPUTS:
%   roc
%   labels_results -- the percentage of dipoles in each label
%_______________________________________________________
% COORDINATES INFORMATION
%   GridLoc: [Nvertices x 3], (x,y,z) positions of the grid of source points.
%   In the case of a surface head model, it corresponds to a copy of the
%   'Vertices' matrix from the cortex surface file.
%
%   https://neuroimage.usc.edu/brainstorm/Tutorials/HeadModel
%
% ATLAS INFORMATION
%   To have cortical and subrotical surfaces in head Model this tutorial
%   should be done: https://neuroimage.usc.edu/brainstorm/Tutorials/DeepAtlas
%   see also: https://neuroimage.usc.edu/brainstorm/Tutorials/ExploreAnatomy

%% Atlas selection
for i = 1:numel(atlas)
    if (strcmp(atlas(i).Name, 'Destrieux'))
        cortex_atlas = atlas(i).Scouts;
    end
end

%% Create labels_results
labels_results = {'Lables','v_mag','v_grad','v_both','mag','grad','both'};
labels_results(2:size({cortex_atlas.Label},2)+1,1) = {cortex_atlas.Label}';



%% Fill labels_results and cluster_mask

labels_results = fill_labels_results(labels_results, 2, v_mag(:,7)', cortex_atlas);
labels_results = fill_labels_results(labels_results, 3, v_grad(:,7)', cortex_atlas);
labels_results = fill_labels_results(labels_results, 4, [v_mag(:,7)', v_grad(:,7)'], cortex_atlas);

labels_results = fill_labels_results(labels_results, 5, mag(:,7)', cortex_atlas);
labels_results = fill_labels_results(labels_results, 6, grad(:,7)', cortex_atlas);
labels_results = fill_labels_results(labels_results, 7, [mag(:,7)', grad(:,7)'], cortex_atlas);

%% Calculate ROC

roc = {'Channels','Number_of_spikes' 'TP', 'FP', 'FN', 'TN','Sensitivity','Specificity'};
roc(2:4,1) = {'grad','mag','both'};
roc(2:4,2) = {length(grad),length(mag),length(grad)+length(mag)};

results = cell2mat(labels_results(2:end,2:end));
    
% Calculate true positive
roc(2,3) = {size(results(results(:,1)>0 & results(:,4)>0),1)};
roc(3,3) = {size(results(results(:,2)>0 & results(:,5)>0),1)};
roc(4,3) = {size(results(results(:,3)>0 & results(:,6)>0),1)};

% Calculate false negative
roc(2,5) = {size(results(results(:,1)>0 & results(:,4)==0),1)};
roc(3,5) = {size(results(results(:,2)>0 & results(:,5)==0),1)};
roc(4,5) = {size(results(results(:,3)>0 & results(:,6)==0),1)};

% Sensitivity
roc(2,7) = {size(results(results(:,1)>0 & results(:,4)>0),1)/size(results(results(:,1)>0),1)};
roc(3,7) = {size(results(results(:,2)>0 & results(:,5)>0),1)/size(results(results(:,2)>0),1)};
roc(4,7) = {size(results(results(:,3)>0 & results(:,6)>0),1)/size(results(results(:,3)>0),1)};

% Calculate false positive
roc(2,4) = {size(results(results(:,1)==0 & results(:,4)>0),1)};
roc(3,4) = {size(results(results(:,2)==0 & results(:,5)>0),1)};
roc(4,4) = {size(results(results(:,3)==0 & results(:,6)>0),1)};

% Calculate true negative
roc(2,6) = {size(results(results(:,1)==0 & results(:,4)==0),1)};
roc(3,6) = {size(results(results(:,2)==0 & results(:,5)==0),1)};
roc(4,6) = {size(results(results(:,3)==0 & results(:,6)==0),1)};

% Specificity

roc(2,8) = {size(results(results(:,1)>0 & results(:,4)==0),1)/size(results(results(:,1)==0),1)};
roc(3,8) = {size(results(results(:,2)>0 & results(:,5)==0),1)/size(results(results(:,2)==0),1)};
roc(4,8) = {size(results(results(:,3)>0 & results(:,6)==0),1)/size(results(results(:,3)==0),1)};

end




%% Find Label by vetex

function label_index = find_label(atlas, vertex)
fun = @(x) sum(atlas(x).Vertices == vertex);
true_label = arrayfun(fun, 1:numel(atlas),'UniformOutput',false);
true_label_vector = cell2mat(true_label)';
label_index = find(true_label_vector);
end

%% Fill labels_results

function labels_results = fill_labels_results(labels_results,...
    column, Vertices, atlas)

results = zeros(size(labels_results,1),1);
for i = 1:length(Vertices)
    
    label_index = find_label(atlas, Vertices(i)); %one vertex could be in two labels
    for ii = label_index
        %!!! ii +1 because labels in labels_results started from 2 !!!
        results(ii+1,1) = results(ii+1,1) + 1;
    end
end
results(:,1) = results(:,1)/sum(results(:,1));
labels_results(2:end,column) = num2cell(results(2:end,1));
end



