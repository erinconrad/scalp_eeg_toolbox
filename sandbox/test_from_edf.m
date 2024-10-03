
%{

notes:
EMU0766_Day02_1_spike_at_79860.21875 nice example


EMU1614_Day05_1_spike_at_62499.75.edf one with probs
EMU0766_Day01_1_spike_at_18757.0625 SEIZURE
%}
clear

save_plot = 0;

%filepath = '/Users/erinconrad/Desktop/research/epilepsy_prediction/Brian Prager results/sub-00001_ses-preimplant001_task-task_run-03_eeg.edf';
%folder_path = '/Users/erinconrad/Library/CloudStorage/Box-Box/EDF_SpikeNet2/';
folder_path = '/Users/erinconrad/Library/CloudStorage/Box-Box/SN12_three_threshold/SN2/edf_threshold_8/';
filename = 'EMU1362_Day02_1_spike_at_2647.6875';
filepath = [folder_path,filename,'.edf'];
times = [0 15];

%% File locs and set path
locations = scalp_toolbox_locs;
% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));


%% Get ieeg data
[values,chLabels,fs] = read_in_edf(filepath);
data.fs = fs;
data.values = values;

%samples = 

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
test_values = bandstop(values,[58 62],fs);

% Demean the channels
values = values - nanmean(values,1);

[bipolar_values,bipolar_labels] = scalp_bipolar(chLabels,values);
data.bipolar_values = bipolar_values;
data.bipolar_labels = bipolar_labels;
plot_scalp_eeg_gain(bipolar_values,fs,bipolar_labels,20);
%plot_scalp_eeg(bipolar_values,fs,bipolar_labels);

if save_plot
    print(gcf,['../../results/',filename,'.png'],'-dpng')
end