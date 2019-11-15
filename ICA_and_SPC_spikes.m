function ICA_and_SPC_spikes(path_cluster_out)
%
% -------------------------------------------------------------------------
% Computation and saving overlap between detections
% -------------------------------------------------------------------------
% ICA_grad
% ICA_mag
% SPC_grad
% SPC_mag
%
% OUTPUTS:
%
% cluster_out_overlap_grad.mat
% cluster_out_overlap_mag.mat
% cluster_out_overlap_grad.csv
% cluster_out_overlap_mag.csv
%
%_______________________________________________________
%

% Load Spyking Circus detected timestamps
SPC_grad = load([path_cluster_out 'SpyCir_based_grad.csv']);
SPC_mag = load([path_cluster_out 'SpyCir_based_mag.csv']);

% Load ICA detected timestamps
ICA_grad = load([path_cluster_out 'ICA_based_grad.csv']);
ICA_mag = load([path_cluster_out 'ICA_based_mag.csv']);

overlap = [];

win = 20; % window around the detection (ms)


for i = 1:length(SPC_grad)
	[nearest, ~] = find_nearest(ICA_grad(:,1), SPC_grad(i,1));
	if (abs(nearest - SPC_grad(i,1))<win)
		overlap = [overlap; SPC_grad(i,:)];
	end
end

csvwrite([path_cluster_out 'overlap_grad.csv'], overlap);
save([path_cluster_out 'overlap_grad.csv'], 'overlap');

overlap = [];
for i = 1:length(SPC_mag)
	[nearest, ~] = find_nearest(ICA_mag(:,1), SPC_mag(i,1));
	if (abs(nearest - SPC_mag(i,1))<win)
		overlap = [overlap; SPC_mag(i,:)];
	end
end
csvwrite([path_cluster_out 'overlap_mag.csv'], overlap);
save([path_cluster_out 'overlap_mag.csv'], 'overlap');

end

function [nearest, index] = find_nearest(stamps, value)
    [~,I] = min(abs(stamps-value));
    nearest = stamps(I);
    index = I;
end



