

clear;
clc;

%% Parameters
just_psd = 0;
psd_ch = 'P3-O1';
psd_times = [7*15 8*15];
duration = 15;
start_time = 0;

%% Paths - edit this!!!
%edf_file = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/temple_test/aaaaacad_s003_t000.edf';
edf_path = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/temple_test/';
spec_file = 'aaaaacad_s003_t000.edf';
edf_file = fullfile(edf_path,[spec_file]);

%% Load and process the EDF file
[values, chLabels, fs] = read_in_edf_clean(edf_file); % Ensure this function is available
chLabels = clean_temple(chLabels);

% Initialize gain
gain = 50;

% Interpolate over NaNs
for j = 1:size(values, 2)
    nan_indices = isnan(values(:, j));
    if any(nan_indices)
        values(nan_indices, j) = nanmean(values(:, j));
    end
end

% Demean the channels
values = values - nanmean(values, 1);

% Make channels with all nans be zero
all_nan_chs = all(isnan(values),1);
values(:,all_nan_chs) = zeros(size(values,1),sum(all_nan_chs));

% High-pass filter
%values = highpass(values, 1, fs);
values = bandstop(values, [58 62], fs);

% Get bipolar montage
[bipolar_values, bipolar_labels] = scalp_bipolar_clean(values, chLabels);

% Also get CAR montage
[car_values, car_labels] = car_montage_clean(values, chLabels);

%% Plotting
% Initially set bipolar montage
montage_type = 'bipolar';
switch montage_type
    case 'bipolar'
        plot_values = bipolar_values;
        plot_labels = bipolar_labels;
    case 'car'
        plot_values = car_values;
        plot_labels = car_labels;
end

if just_psd
    figure
    set(gcf,'position',[1 1 1400 500])
    tiledlayout(2,1)

    psd_idx =  [max([round(psd_times(1)*fs),1]):round(psd_times(2)*fs)];
    temp_values = bipolar_values(psd_idx,strcmp(bipolar_labels,psd_ch));

    nexttile
    plot(linspace(0,psd_times(2)-psd_times(1),length(temp_values)),temp_values)
    hold on
    for k = 1:15
        plot([k k],get(gca,'ylim'),'r--')
    end
    xlabel('seconds')
    ylabel('uV)')
    title(sprintf('%s %s',spec_file,psd_ch))

    nexttile
    x = temp_values;
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    %psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx = log(abs(xdft).^2);
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/length(x):fs/2;
    plot(freq,psdx);
    xlim([0 70])
    xlabel('Hz')
    ylabel('log(uV^2)')
    
    print(gcf,[edf_path,strrep(spec_file,'.',''),'_',num2str(psd_times(1)),'-',num2str(psd_times(2))],'-dpng')


else

    % Create a figure for displaying the EEG data
    figure;
    hFig = gcf;
    set(hFig, 'Position', [10 10 1400 1000]);
    tiledlayout(1, 1, 'Padding', 'compact');
    nexttile;
    
    % Initialize user_input
    user_input = '';  % Initialize as empty
    
    % Plot EEG data with prompt and selection texts
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);
    
    fprintf('\nDisplaying start time %1.1f s\n', start_time);
    
    %% Store Handles and Data in appdata
    setappdata(hFig, 'gain', gain);
    setappdata(hFig, 'bipolar_values', bipolar_values);
    setappdata(hFig, 'bipolar_labels', bipolar_labels);
    setappdata(hFig, 'car_values', car_values);
    setappdata(hFig, 'car_labels', car_labels);
    setappdata(hFig, 'fs', fs);
    setappdata(hFig, 'montage_type', montage_type);
    setappdata(hFig, 'start_time', start_time);
    setappdata(hFig, 'duration', duration);
    
    %% Handle User Input
    % Initialize variables and store them in appdata
    arrow_key_input = '';
    setappdata(hFig, 'arrow_key_input', arrow_key_input);
    
    % Set up the figure to detect key presses using WindowKeyPressFcn
    set(hFig, 'WindowKeyPressFcn', @keyPressHandler);
    
    % Wait until the user has confirmed their response
    uiwait(hFig);  % Wait until uiresume is called
    
    %% Retrieve User Input and Other Data
    % Retrieve the stored arrow_key_input and user_input from appdata
    arrow_key_input = getappdata(hFig, 'arrow_key_input');
    gain = getappdata(hFig, 'gain');
    montage_type = getappdata(hFig, 'montage_type');
    
    %% Handle End of File and Write Responses
    % Now you can safely close the figure
    if ishandle(hFig)
        close(hFig);
    end

end

%% Function Definitions
function new_labels = clean_temple(chLabels)

new_labels = chLabels;
new_labels = strrep(new_labels,'EEG','');
new_labels = strrep(new_labels,'_REF','');
new_labels = strrep(new_labels,'FP','Fp');
new_labels = strrep(new_labels,'Z','z');

end


function [out, chLabels, fs] = read_in_edf_clean(filepath)

data = edfread(filepath);
info = edfinfo(filepath);
fs = info.NumSamples(1);

chs = info.SignalLabels;
nchs = length(chs);
chs = replace(chs, "-", "_");
chs = replace(chs, " ", "");

nsamples = size(data,1) * fs;

out = nan(nsamples,nchs);

for i = 1:nchs
    if iscell(data.(chs(i)))
        out(:,i) = cell2mat(data.(chs(i)));
    else
        continue
    end

end

chLabels = cellstr(chs);

end



function keyPressHandler(hObject, event)
    % Nested key press handler function

    % Retrieve data stored in appdata
    arrow_key_input = getappdata(hObject, 'arrow_key_input');
    gain = getappdata(hObject, 'gain');
    bipolar_values = getappdata(hObject, 'bipolar_values');
    bipolar_labels = getappdata(hObject, 'bipolar_labels');
    car_values = getappdata(hObject, 'car_values');
    car_labels = getappdata(hObject, 'car_labels');
    fs = getappdata(hObject, 'fs');
    montage_type = getappdata(hObject, 'montage_type');
    start_time = getappdata(hObject, 'start_time');
    duration = getappdata(hObject, 'duration');

    % Handle key presses
    if strcmp(event.Key, 'uparrow')
        disp('Up arrow pressed');
        arrow_key_input = 'up';
        gain = gain - 10;
        setappdata(hObject, 'gain', gain);
        hold off;
        % Re-plot with updated gain
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end

        % Re-plot EEG with updated gain and existing user_input
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);

    elseif strcmp(event.Key, 'downarrow')
        disp('Down arrow pressed');
        arrow_key_input = 'down';
        gain = gain + 10;
        setappdata(hObject, 'gain', gain);
        hold off;
        % Re-plot with updated gain
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        % Re-plot EEG with updated gain and existing user_input
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);

     
    elseif strcmp(event.Key, 'leftarrow')
        disp('Left arrow pressed')
        arrow_key_input = 'left';
        start_time = max([start_time - duration,0]);
        setappdata(hObject,'start_time',start_time);
        hold off;
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end

        % Re-plot EEG with updated gain and existing user_input
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);

    elseif strcmp(event.Key, 'rightarrow')
        disp('Right arrow pressed')
        arrow_key_input = 'right';
        start_time = min([start_time + duration,size(bipolar_values,1)/fs]);
        setappdata(hObject,'start_time',start_time);
        hold off;
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end

        % Re-plot EEG with updated gain and existing user_input
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);

    elseif strcmp(event.Key, 'b') || strcmp(event.Key, 'c')
        disp(['Montage change requested: ', event.Key]);
        if strcmp(event.Key, 'b')
            montage_type = 'bipolar';
        elseif strcmp(event.Key, 'c')
            montage_type = 'car';
        end
        setappdata(hObject, 'montage_type', montage_type);
        hold off;
        % Re-plot with updated montage
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain, start_time,duration);

      
    
    else
        disp(['Other key pressed: ', event.Key]);
    end

    % Store the updated data back into appdata
    setappdata(hObject, 'arrow_key_input', arrow_key_input);
end

function [bipolar_values, bipolar_labels] = scalp_bipolar_clean(values, chLabels)
    % This function creates a bipolar montage
    
    montage = {'Fp1','F7';...
        'F7','T3';...
        'T3','T5';...
        'T5','O1';...
        '','';...
        'Fp2','F8';...
        'F8','T4';...
        'T4','T6';...
        'T6','O2';...
        '','';...
        'Fp1','F3';...
        'F3','C3';...
        'C3','P3';...
        'P3','O1';...
        '','';...
        'Fp2','F4';...
        'F4','C4';...
        'C4','P4';...
        'P4','O2';...
        '','';...
        'Fz','Cz';...
        '','';...
        '','';...
        'EKG1',''};
    
    bipolar_labels = cellfun(@(x,y) sprintf('%s-%s',x,y), montage(:,1), montage(:,2), 'UniformOutput', false);
    
    nsamples = size(values,1);
    nbi_channels = size(montage,1);
    bipolar_values = nan(nsamples, nbi_channels);
    
    % Convert ECGL and ECGR to EKG1 and EKG2
    if sum(strcmp(chLabels, 'ECGL')) ~= 0
        chLabels(strcmp(chLabels, 'ECGL')) = {'EKG1'};
    end
    if sum(strcmp(chLabels, 'ECGR')) ~= 0
        chLabels(strcmp(chLabels, 'ECGR')) = {'EKG2'};
    end
    
    for ib = 1:nbi_channels
        
        curr_mon = montage(ib,:);
        
        % Leave the empty rows as NaNs (just to designate blank space)
        if strcmp(curr_mon{1}, '') && strcmp(curr_mon{2}, '')
            continue;
        end
        
        % Get the channel labels that match
        lab1 = strcmp(chLabels, curr_mon{1});
        lab2 = strcmp(chLabels, curr_mon{2});
    
        if ~(sum(lab1) == 1 && sum(lab2) == 1)
            if contains(curr_mon{1}, 'EKG')
                if sum(lab1) == 0
                    continue;
                end
            else
                error('Missing expected channel: %s or %s', curr_mon{1}, curr_mon{2});
            end
        end
    
        if sum(lab2) == 0 && strcmp(chLabels(lab1), 'EKG1')
            bipolar_values(:, ib) = values(:, lab1);
            continue;
        end
        
        % Get difference in EEG signal
        val1 = values(:, lab1);
        val2 = values(:, lab2);
        bipolar_values(:, ib) = val1 - val2;
        
    end
end

function [car_values, car_labels] = car_montage_clean(values, chLabels)
    % This function creates a CAR (Common Average Reference) montage
    
    car_labels = chLabels;
    brain_channels = {'Fp1',...
        'F7',...
        'T3',...
        'T5',...
        'O1',...
        'Fp2',...
        'F8',...
        'T4',...
        'T6',...
        'O2',...
        'F3',...
        'C3',...
        'P3',...
        'F4',...
        'C4',...
        'P4',...
        'Fz',...
        'Cz'};
    
    % Take average of the brain channels
    is_brain = ismember(chLabels, brain_channels);
    avg = nanmean(values(:, is_brain), 2);
    
    % Do CAR (just average the non-skip channels)
    car_values = values; % Initialize
    car_values(:, is_brain) = values(:, is_brain) - repmat(avg, 1, sum(is_brain));
    
    for i = 1:length(chLabels)
        if is_brain(i) == 1
            car_labels{i} = [chLabels{i}, '-CAR'];
        end
    end
end

function plot_scalp_eeg_clean(values, fs, labels, gain,start_time,duration)
    % This function plots the EEG data in a scalp montage
    
    added_offset = gain;
    
    dur = duration;
    nchs = length(labels);

    % restrict values
    start_idx = max([1,round(start_time*fs)]);
    end_idx = min([size(values,1),round((start_time+dur)*fs)]);
    values = values(start_idx:end_idx,:);
    
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
