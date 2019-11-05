function roc = roc_time(visual, grad, mag)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% Sensitivity and Specificity for comparison with visual detections
%  
% -------------------------------------------------------------------------
%
% USAGE:
% roc = roc_curve(visual, grad, mag)
%
% PARAMETERS:
% win - the window around the visual spike (+- win) in which a detection will be marked as True
%
% INPUT:
% Read data as matrices with two clumns
% Time should be in ms
% visual - one column (timestamps)
% grad - two columns (timestamps, clusters)
% mag - two columns (timestamps, clusters)
%
% OUTPUT:
% roc - table with columns: 'Channels','Number_of_spikes' 'TP', 'FP', 'FN', 'TN','Sensitivity','Specificity'
%                  rows: 'grad','mag','both'
% -------------------------------------------------------------------------


    roc = {'Channels','Number_of_spikes' 'TP', 'FP', 'FN', 'TN','Sensitivity','Specificity'};
    roc(2:4,1) = {'grad','mag','both'};
    roc(2:4,2) = {length(grad),length(mag),length(grad)+length(mag)};
    
    win = 20; % window around visual spike (ms)
    mag(:,3) = 0;
    grad(:,3) = 0;

%% 1 Sensitivity
    % columns in visual array:
    % 1 - time of visual spikes
    % 2 - time of nearest grad spike
    % 3 - cluster of nearest grad spike
    % 4 - absolute difference between visual spike and grad
    % 5 - time of nearest mag spike
    % 6 - cluster of nearest mag spike
    % 7 - absolute difference between visual spike and mag
    % 8 - if 0 - false positive, if 1 - true positive (both mag and grad)
    for i = 1:size(visual,1)
        [visual(i,2), ind] = find_nearest(grad(:,1), visual(i,1));
        visual(i,3) = grad(ind,2);
        visual(i,4) = abs(visual(i,1) - grad(ind,1));
        if (visual(i,4)<win) % threshold in ms
            grad(ind,3) = true; % true positive spike and cluster
            visual(i,8) = true;
        end

        [visual(i,5), ind] = find_nearest(mag(:,1), visual(i,1));
        visual(i,6) = mag(ind,2);
        visual(i,7) = abs(visual(i,1) - mag(ind,1));
        if (visual(i,7)<win) % threshold in ms
            mag(ind,3) = true; % true positive spike
            visual(i,8) = true;
        end
    end

    % Adding clusters to true positive in which TP spikes were found
    for cluster = 1:max(grad(:,2))
        ind_cluster = grad(:,2) == cluster;
        if sum(grad(ind_cluster,3))>0
            grad(ind_cluster,3) = true;
        end
    end
    for cluster = 1:max(mag(:,2))
        ind_cluster = mag(:,2) == cluster;
        if sum(mag(ind_cluster,3))>0
            mag(ind_cluster,3) = true;
        end
    end

    % Calculate true positive
    roc(2,3) = {length(visual(visual(:,4)<win,1))};
    roc(3,3) = {length(visual(visual(:,7)<win,1))};
    roc(4,3) = {length(visual((visual(:,4)<win)|(visual(:,7)<win),1))};

    % Calculate false negative
    roc(2,5) = {length(visual(visual(:,4)>=win,1))};
    roc(3,5) = {length(visual(visual(:,7)>=win,1))};
    roc(4,5) = {length(visual(visual(:,4)>=win,1))+length(visual(visual(:,7)>=win,1))};

    % Calculate sensitivity
    roc(2,7) = {length(visual(visual(:,4)<win,1))/length(visual)};
    roc(3,7) = {length(visual(visual(:,7)<win,1))/length(visual)};
    roc(4,7) = {length(visual((visual(:,4)<win)|(visual(:,7)<win),1))/length(visual)};

%% 2 Specificity
    % Calculate false positive

    roc(2,4) = {length(grad(grad(:,3)==false,1))};
    roc(3,4) = {length(mag(mag(:,3)==false,1))};
    roc(4,4) = {length(grad(grad(:,3)==false,1))+length(mag(mag(:,3)==false,1))};

    % Calculate true negative
    %TN_points(:,1) - time points without visual spike after them in interval 100ms
    %TN_points(:,2) - true negative grad
    %TN_points(:,3) - true negative mag

    TN_points = zeros(length(visual(:,1)),3);
    row = 1;
    while length(TN_points(TN_points(:,1) ~= 0,1)) < length(visual(:,1))
%         rt = randi([1 max(visual(:,1))-100],1,1); %random time point
        rt = randi([1 60000-100],1,1); %random time point

        if ((length(visual((visual(:,1)>rt)&(visual(:,1)<rt+100),1)))==0)
            TN_points(row,1) = rt;

            if ((length(grad((grad(:,1)>rt)&(grad(:,1)<rt+100),1)))==0)
                TN_points(row,2) = 1;
            end
            if ((length(mag((mag(:,1)>rt)&(mag(:,1)<rt+100),1)))==0)
                TN_points(row,3) = 1;
            end

            row = row + 1;
        end
    end

    roc(2,6) = {length(TN_points(TN_points(:,2)==true,1))};
    roc(3,6) = {length(TN_points(TN_points(:,3)==true,1))};
    roc(4,6) = {length(TN_points(((TN_points(:,2)==true)&(TN_points(:,3)==true)),1))};

    % Specificity
    TN = roc(2,6);
    FP = roc(2,4);
    roc(2,8) = {TN{:}/(TN{:}+FP{:})};
    TN = roc(3,6);
    FP = roc(3,4);
    roc(3,8) = {TN{:}/(TN{:}+FP{:})};
    TN = roc(4,6);
    FP = roc(4,4);
    roc(4,8) = {TN{:}/(TN{:}+FP{:})};


end
%% Find nearest spike function
    % time - column of timestamps
    % value - one time point

function [nearest, index] = find_nearest(stamps, value)
    [~,I] = min(abs(stamps-value));
    nearest = stamps(I);
    index = I;
end

