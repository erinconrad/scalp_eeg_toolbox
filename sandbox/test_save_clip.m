function test_save_clip

file_name = 'test_scalp_withcm';
start_time = 71485.33;
end_time = start_time + 15;

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

[bipolar_values,bipolar_labels] = scalp_bipolar(chLabels,values);
data.bipolar_values = bipolar_values;
data.bipolar_labels = bipolar_labels;
plot_scalp_eeg(bipolar_values,fs,bipolar_labels);

%% Save the clip
save([locations.main_folder,'data/eeg_clips/',file_name],'data');

end