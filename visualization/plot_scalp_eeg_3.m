function plot_scalp_eeg_3(values, fs, labels, gain)
    % This function plots the EEG data in a scalp montage
    
    added_offset = gain;
    
    dur = size(values,1)/fs;
    nchs = length(labels);
    
    offset = 0;
    
    ch_offsets = zeros(nchs,1);
    ch_bl = zeros(nchs,1);
    
    for ich = 1:nchs
        plot(linspace(0, dur, size(values,1)), values(:, ich) - offset, 'k');
        hold on;
        ch_offsets(ich) = offset;
        ch_bl(ich) = -offset + nanmedian(values(:, ich));
        offset = offset + added_offset;
        text(dur + 0.05, ch_bl(ich), sprintf('%s', labels{ich}), 'FontSize', 15);
    end
    plot([dur/2 dur/2], get(gca, 'ylim'), 'r--');
    
    yticklabels([]);
    xticks(1:floor(dur));
    ylim([-offset - added_offset, added_offset]);
    xlabel('Time (seconds)');
    set(gca, 'FontSize', 15);
    
    %% Add second markers
    for i = 1:floor(dur)
        plot([i, i], get(gca, 'ylim'), 'k--');
    end
end