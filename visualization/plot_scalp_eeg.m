function plot_scalp_eeg(values,fs,labels)

added_offset = 20;
figure

set(gcf,'position',[10 10 1400 1000])
tiledlayout(1,1,'padding','compact')
nexttile
dur = size(values,1)/fs;
nchs = length(labels);

offset = 0;

ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for ich = 1:nchs
    plot(linspace(0,dur,size(values,1)),values(:,ich) - offset,'k');
    hold on
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    offset = offset+added_offset;
    text(dur+0.05,ch_bl(ich),sprintf('%s',labels{ich}),'fontsize',15)
    
    if ich < nchs
        if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
            offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
        end
        
    end
    
end

yticklabels([])
xticks(1:floor(dur))
ylim([-offset-added_offset added_offset])
xlabel('Time (seconds)')
set(gca,'fontsize',15)

%% Add second markers
for i = 1:dur
    plot([i,i],get(gca,'ylim'),'k--');
    
end


end