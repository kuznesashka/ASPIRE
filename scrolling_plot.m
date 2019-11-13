 function scrolling_plot(x, Y, num_val, events, nw, kurtosis, channels, channel_idx, W)
% -------------------------------------------------------
% Visualization of the ICA timeseries sorted in the
% order of the kurtosis decrease, spikes visualization
% on the second graph
% -------------------------------------------------------
% FORMAT:
%   scrolling_plot(x, Y, num_val, events, nw, kurtosis)
% INPUTS:
%   x -- [1 x T] values on the x axis
%   Y -- [ncomps x T] matrix of values on the y axis
%   num_val -- number of values shown at one time, step between the two
%   figures is 0.01*T
%   events -- [1 x T] indicator of the event
%   nw -- number of components shown on the one plot, ncomps::nw
%
% NOTE:
% 
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

% for i = 1:length(channel_idx)
%     chan_loc(i,:) = mean(channel.Channel(channel_idx(i)).Loc, 2)';
% end
% [x_coord,y_coord] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), ...
%     chan_loc(:,3),'2dlayout');


% Create the first figure
f = figure;
% pick the first nw components
y = Y((1:nw),:);
k = kurtosis(1:nw);
for i = [1:nw]
    subplot((nw+1),1,i)
    plot(x(1:(1+num_val)),y(i,1:(1+num_val)))
    xlim([x(1) x(1+num_val)])
    ylim([min(y(i,:)) max(y(i,:))])
    title(['kurtosis = ', num2str(k(i))])
end
% l = 1;
% for i = [2:2:nw*2]
%     subplot((nw+1),2,i)
%     plot_topo(x_coord, y_coord, Q(:,l));
%     l = l+1;
% end
subplot(nw+1,1,(nw+1))
plot(x(1:(1+num_val)), events(1:(1+num_val)),'LineWidth', 1.5)
xlim([x(1) x(1+num_val)])
ylim([0 1])


% Create slide slide the time axis
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',size(y,2)-num_val-mod((size(y,2)-num_val), 100),...
    'Value',1,...
    'Position', [400 5 120 20],...
    'SliderStep', [0.01 1], ... % step = 0.01*T
    'Callback', @slidedata);

% Create pop-up menu to change components
ncomp = size(Y,1)/nw;
p = 1;
for i = 1:ncomp
    names{p} = [num2str((i-1)*nw+1),'-', num2str(i*nw), ' component'];
    p = p+1;
end
popup = uicontrol('Style', 'popup',...
       'String', {char(names(1,:))},...
       'Position', [100 5 120 20],...
       'Callback', @setcomponent);

% Add a text uicontrol to label the slider.
%     txt = uicontrol('Style','text',...
%         'Position',[400 45 120 20],...
%         'String','Vertical Exaggeration');

% Make figure visble after adding all components
f.Visible = 'on';

function slidedata(source,event)
    val = int32(source.Value);
    for i = 1:nw
        subplot(nw+1,1,i)
        plot(x(val:(val+num_val)),y(i,val:(val+num_val)))
        xlim([x(val) x(val+num_val)])
        ylim([min(y(i,:)) max(y(i,:))])
        title(['kurtosis = ', num2str(k(i))])
    end
    subplot(nw+1,1,(nw+1))
    plot(x(val:(val+num_val)), events(val:(val+num_val)), 'LineWidth', 1.5)
    xlim([x(val) x(val+num_val)])
    ylim([0 1])
 end

function setcomponent(source,event)
    valcomp = source.Value;
    y = Y((((valcomp-1)*nw+1):valcomp*nw),:);
    k = kurtosis((((valcomp-1)*nw+1):valcomp*nw));
   
    val = int32(get(sld,'value'));
    for i = 1:nw
        subplot(nw+1,1,i)
        plot(x(val:(val+num_val)),y(i,val:(val+num_val)))
        xlim([x(val) x(val+num_val)])
        ylim([min(y(i,:)) max(y(i,:))])
        title(['kurtosis = ', num2str(k(i))])
    end
    subplot(nw+1,1,(nw+1))
    plot(x(val:(val+num_val)), events(val:(val+num_val)), 'LineWidth', 1.5)
    xlim([x(val) x(val+num_val)])
    ylim([0 1])        
end

end