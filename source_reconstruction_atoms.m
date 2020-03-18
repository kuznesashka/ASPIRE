function VE = source_reconstruction_atoms(Data, G3, channel_idx, dip_ind, N)

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
    VE 		= zeros(length(dip_ind), length(Data)); 
    chunk_len = length(Data)/N;
    for chunk = 1:length(VE)/chunk_len
        chunk_begin = 1 + (chunk-1)*chunk_len;
        chunk_end   = chunk_len + (chunk-1)*chunk_len;
        chunk_data  = Data(:, chunk_begin:chunk_end);
        % tic
        for i = 1:length(dip_ind) % all atoms locations
            g = G2(:,(dip_ind(i)*2-1):dip_ind(i)*2); %[sensors x spikes_locationx2]
            [U S ~] = svd(chunk_data); % U:[204 x 204], S:[204xtime]
            h = cumsum(diag(S)/sum(diag(S)));
            n = min(find(h>=0.95));
            [u s v] = svd(U(:,1:n)'*g); % project only main components, SVD them again: V is the oreintation of the source
            g_fixed = g*v(1,:)'; % fixed orientation forward model
            VE(i, chunk_begin:chunk_end) = chunk_data'*g_fixed; % estimated source signal
        end
        % toc
    end
end