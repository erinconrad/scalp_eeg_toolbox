%% example
T = readtable('../../results/EMU1614_Day05_1_time_probabilities_Sep2924.csv');

folder_path = '/Users/erinconrad/Library/CloudStorage/Box-Box/EDF_SpikeNet2/';
%filepath = [folder_path,'EMU1614_Day05_1_spike_at_62499.75.edf'];
filepath = ['../../results/EMU1614_Day05_1_62492_to_62507.edf'];
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

figure
tiledlayout(2,1,'TileSpacing','compact')
nexttile
plot(linspace(0,15,size(bipolar_values,1)),bipolar_values(:,7))
hold on
plot([7.5 7.5],ylim,'k--')
ylabel('Amplitude on F8-T4 (uV)')
xticklabels([])
nexttile
plot(linspace(0,15,length(T.PredictedProbability)),T.PredictedProbability)
hold on
plot([7.5 7.5],ylim,'k--')
ylabel('SpikeNet2 probability')
xlabel('Time (s)')
print(gcf,'../../results/spikeNet_example','-dpng')