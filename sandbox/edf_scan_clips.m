
%{

notes:
EMU0766_Day02_1_spike_at_79860.21875 nice example


EMU1614_Day05_1_spike_at_62499.75.edf one with probs
EMU0766_Day01_1_spike_at_18757.0625 SEIZURE
%}
clear

save_plot = 0;


%folder_path = '../../data/eeg_clips/';
%filename = 'EMU1361_Day01_1_4042.96875_to_4429.6875';
folder_path = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/juri_code/';
filename = 'EMU1371_Day02_1_5006_to_5491.edf';

filepath = [folder_path,filename];
clip_dur = 15; % how many seconds to display at a time
%times = [350-7.5 350+7.5];

%% File locs and set path
locations = scalp_toolbox_locs;
% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Get eeg data
[values,chLabels,fs] = read_in_edf(filepath);


data.fs = fs;
data.values = values;
eeg_times = linspace(0,size(values,1)/fs,size(values,1));

if sum(~isnan(values),'all') == 0
    error('Data requested is all nans!')
end

% Interpolate over nans
for i = 1:size(values,2)
    values(isnan(values(:,i)),i) = nanmean(values(:,i));
end

% High pass filter
old_values = values;
values = highpass(old_values,0.5,fs);
values = bandstop(values,[58 62],fs);

% Demean the channels
values = values - nanmean(values,1);

% Loop over times to just show requested duration
ntimes = ceil(eeg_times(end)/clip_dur);
nsamples_per_clip = round(clip_dur * fs);
for it = 1:ntimes
    
    curr_samples = nsamples_per_clip*(it-1)+1:min(nsamples_per_clip*it,size(values,1));
    curr_times = eeg_times(curr_samples);
    curr_values =values(curr_samples,:);
    [bipolar_values,bipolar_labels] = scalp_bipolar(chLabels,curr_values);

    figure
    set(gcf,'Position',[1 1 1400 1000])
    plot_scalp_eeg_gain(bipolar_values,fs,bipolar_labels,50);
    title(sprintf('Clip %d',it))

    fprintf('\nShowing clip %d of %d. Press any button to move to next one.\n',...
        it,ntimes)
    pause
    close gcf

end




if save_plot
    print(gcf,['../../results/',filename,'.png'],'-dpng')
end