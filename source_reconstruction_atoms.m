function VE = source_reconstruction_atoms(Data, G3, channel_idx, dip_ind)

% -------------------------------------------------------------------------
% Source timeseries (dipole-fitting)
% -------------------------------------------------------------------------
% INPUTS:
%   Data -- brainstorm structure with artifact corrected maxfiltered MEG data
%   G3 -- brainstorm structure with forward operator
%   channel_type -- channels used ('grad' or 'mag')
%   
% OUTPUTS:
%   VE -- estimated source signal
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

    % 2D forward operator
    [G2, ~] = G3toG2(G3, channel_idx);
    VE 		= zeros(length(Date), length(dip_ind)) 

    for i = 1:length(dip_ind) % all atoms locations
        g = G2(:,(dip_ind(i)*2-1):dip_ind(i)*2); %[sensors x spikes_locationx2]
        [U S ~] = svd(Data); % U:[204 x 204], S:[204xtime]
        h = cumsum(diag(S)/sum(diag(S)));
        n = min(find(h>=0.95));
        [u s v] = svd(U(:,1:n)'*g); % project only main components, SVD them again: V is the oreintation of the source
        g_fixed = g*v(1,:)'; % fixed orientation forward model
        VE(:,i) = Data'*g_fixed; % estimated source signal
    end    
end