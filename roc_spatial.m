function roc = roc_spatial(visual_grad, visual_mag, grad, mag)
%
% -------------------------------------------------------------------------
% Spatial ROC curve
% -------------------------------------------------------------------------
% INPUTS:
%   visual_grad -- 3 coordinats of manual dipoles based on grad
%   visual_mag -- 3 coordinats of manual dipoles based on mag
%   grad -- 3 coodinats of dipoles based on grad
%   mag -- 3 coodinats of dipoles based on grad
% 
% OUTPUTS:
%   roc -- tble with results
%
% PARAMETER:
%   distance_threshold -- distance between the manual dipole and the nearest our dipole
%_______________________________________________________


    roc = {'Channels','Number_of_spikes' 'TP', 'FP', 'FN', 'TN','Sensitivity','Specificity'};
    roc(2:4,1) = {'grad','mag','both'};
    roc(2:4,2) = {length(grad),length(mag),length(grad)+length(mag)};
    % distance between the manual dipole and the nearest our dipole
    distance_threshold = 0.01;

%% 1 Sensitivity

    for i = 1:size(visual_grad,1)
        %grad
        [minDistance, index] = nearest_manual_dipole(grad(:,4:6), visual_grad(i,4:6));   
        if (minDistance<distance_threshold) % threshold in m(?)
            grad(index,7) = true; % true positive dipole
            visual_grad(i,7) = minDistance;
            visual_grad(i,8) = true;
        else
            grad(index,7) = false; % true positive dipole
            visual_grad(i,7) = minDistance;
            visual_grad(i,8) = false;
        end
    end
        %mag
            for i = 1:size(visual_mag,1)

        [minDistance, index] = nearest_manual_dipole(mag(:,4:6), visual_mag(i,4:6));
        if (minDistance<distance_threshold) % threshold in m(?)
            mag(index,7) = true; % true positive dipole
            visual_mag(i,7) = minDistance;
            visual_mag(i,8) = true;
        else
            mag(index,7) = false; % true positive dipole
            visual_mag(i,7) = minDistance;
            visual_mag(i,8) = false;
            
        end
    end

    NN = size(visual_grad,1);
    
    % Calculate true positive
    roc(2,3) = {sum(visual_grad(:,8))};
    roc(3,3) = {sum(visual_mag(:,8))}; 
    roc(4,3) = {sum(visual_grad(:,8) | visual_grad(:,8))};
    
    
        NN_grad = size(visual_grad,1);
        NN_mag = size(visual_mag,1);


    % Calculate false negative
    roc(2,5) = {NN_grad - sum(visual_grad(:,8))};
    roc(3,5) = {NN_mag - sum(visual_mag(:,8))};
    roc(4,5) = {NN_grad+NN_mag - (sum(visual_grad(:,8)) + sum( visual_mag(:,8)))};

    % Calculate sensitivity
    roc(2,7) = {sum(visual_grad(:,8))/NN_grad};
    roc(3,7) = {sum(visual_mag(:,8))/NN_mag};
    roc(4,7) = {(sum(visual_grad(:,8)) + sum( visual_mag(:,8)))/(NN_grad+NN_mag)};

%% 2 Specificity


end


%% Find nearest spike function
    

function [minDistance, index] = nearest_manual_dipole(array, vector)
    distances = pdist2(vector, array);
    [minDistance, index] = min(distances);
end




