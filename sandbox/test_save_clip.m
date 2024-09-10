
%{
file_name = 'test_scalp_withcm';
start_time = 71485.33;
%}

file_name = 'HUP215_phaseII_D03';%'EMU0779_Day01_1';
start_time = 418705.80;%9567.2;
end_time = 418705.80+15;%9567.2+30;

%% File locs and set path
locations = scalp_toolbox_locs;
% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
% ieeg stuff
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Get ieeg data
data = download_ieeg_data(file_name,login_name,pwfile,[start_time end_time],1); % 1 means get lots of data
values = data.values;
chLabels = data.chLabels(:,1);
fs = data.fs;
data.start_time = start_time;
data.end_time = end_time;
data.file_name = file_name;

if sum(~isnan(values),'all') == 0
    error('Data requested is all nans!')
end

% Interpolate over nans
%for i = 1:size(values,2)
%    values(isnan(values(:,i)),i) = nanmean(values(:,i));
%end

% High pass filter
%old_values = values;
%values = highpass(old_values,0.5,fs);
%test_values = bandstop(values,[58 62],fs);

% Demean the channels
%values = values - nanmean(values,1);

[bipolar_values,bipolar_labels] = scalp_bipolar(chLabels,values);
data.bipolar_values = bipolar_values;
data.bipolar_labels = bipolar_labels;
plot_scalp_eeg(bipolar_values,fs,bipolar_labels);

%% Show PSD in P4-02
%{
ch = find(strcmp(bipolar_labels,'Fz-Cz'));
curr_vals = bipolar_values(:,ch);
x = curr_vals;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")
%}
%% Save the clip
%save([locations.main_folder,'data/eeg_clips/',file_name],'data');

