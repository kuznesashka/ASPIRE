Ncomp = size(picked_comp_top,2);
figure,


for comp= 1:Ncomp
    subplot(Ncomp,3,1+3*(comp-1))
        vector = picked_comp_top(:,comp)
        plot_topography_MEG(vector, 'grad', labels)
end

subplot(Ncomp,3,sort([(1:Ncomp)*3-1 (1:Ncomp)*3]))
for j = 1:Ncomp
    plot(Data.Time, picked_components(j,:)-50*j,'k')
    hold on
    sp = round(spike_ind(find(component_indicatior(:,j))));
    scatter(Data.Time(sp), picked_components(j,sp)-50*j,'ro')
end
 
addScrollbar(gca,15)