function [Valmax, Indmax, Sources] = RAP_MUSIC_scan_atoms(spike, ...
    thresh, G3, channel_idx)
% -------------------------------------------------------
% RAP-MUSIC scan
% -------------------------------------------------------
% INPUTS:
%   spike -- Nch x T MEG data with spike
%   Gain -- Nch x Nsources*3 forward operator
%   G2 -- Nch x Nsources*2 forward operator in tangential plane
%   thresh -- minimal considered correlation value
%
% OUTPUT:
%   Valmax -- subcorr values for all found dipoles higher than threshold
%   IndMax -- location of this sources
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

% 2D forward operator
[G2, ~] = G3toG2(G3, channel_idx);
Gain = G3.Gain(channel_idx,:)

[Ns, Nsrc2] = size(G2);
Nsrc = Nsrc2/2;

Valmax = [];
Indmax = [];
Sources = [];
[U,S,V] = svd(spike);
h = cumsum(diag(S)/sum(diag(S)));
n = find(h>=0.95);
corr = MUSIC_scan(G2, U(:,1:n(1)));
[valmax, indmax] = max(corr);
source_ts = source_reconstruction_atom(spike, G2, indmax, channel_idx);
Sources = [Sources source_ts];

while valmax > thresh
    Valmax = [Valmax, valmax];
    Indmax = [Indmax, indmax];
    
    A = Gain(:,(indmax*3-2):indmax*3);
    P = eye(Ns, Ns)-A*inv(A'*A)*A';
    spike_proj = P*spike;
    G_proj = P*Gain;
    Gain = G_proj;
    
    G2 = zeros(Ns,2*Nsrc);
    range = 1:2;
    for i = 1:Nsrc
        g = [G_proj(:,1+3*(i-1)) G_proj(:,2+3*(i-1)) G_proj(:,3+3*(i-1))];
        [u sv v] = svd(g);
        gt = g*v(:,1:2);
        G2(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
        range = range + 2;
    end

    [U,S,V] = svd(spike_proj);
    h = cumsum(diag(S)/sum(diag(S)));
    n = find(h>=0.95);
    corr = MUSIC_scan(G2, U(:,1:n(1)));
    [valmax, indmax] = max(corr);
    source_ts = source_reconstruction_atom(spike_proj, G2, indmax, channel_idx);
    Sources = [Sources source_ts];
end

end


function source_ts = source_reconstruction_atom(Data, G2, dip_ind, channel_idx)
    g = G2(:,(dip_ind*2-1):dip_ind*2); %[sensors x spikes_locationx2]
    [U S ~] = svd(Data); % U:[204 x 204], S:[204xtime]
    h = cumsum(diag(S)/sum(diag(S)));
    n = min(find(h>=0.95));
    [u s v] = svd(U(:,1:n)'*g); % project only main components, SVD them again: V is the oreintation of the source
    g_fixed = g*v(1,:)'; % fixed orientation forward model
    source_ts = Data'*g_fixed; % estimated source signal

end



