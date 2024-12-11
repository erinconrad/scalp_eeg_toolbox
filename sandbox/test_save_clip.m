
%{
file_name = 'test_scalp_withcm';
start_time = 71485.33;
%}

file_name = 'EMU1371_Day02_1';%'EMU0779_Day01_1';
start_time = 5006.53;%9567.2;
end_time = 5491.88;%9567.2+30;

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
%data.start_time = start_time;
%data.end_time = end_time;
%data.file_name = file_name;

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

%{
[bipolar_values,bipolar_labels] = scalp_bipolar(chLabels,values);
data.bipolar_values = bipolar_values;
data.bipolar_labels = bipolar_labels;
plot_scalp_eeg(bipolar_values,fs,bipolar_labels);
%}


%% Save the clip
%save([locations.main_folder,'data/eeg_clips/',file_name],'data');
% Prep header
nchs = length(chLabels);
nsamples = size(values,1);

hdr = edfheader("EDF+");
hdr.NumDataRecords = 1;
hdr.DataRecordDuration = seconds(nsamples/fs);
hdr.NumSignals = nchs;
hdr.SignalLabels = chLabels;
hdr.PhysicalDimensions = repelem("uV",nchs);
hdr.PhysicalMin = min(values);
hdr.PhysicalMax = max(values);

hdr.PhysicalMin(isnan(hdr.PhysicalMin)) = -1e3;
hdr.PhysicalMax(isnan(hdr.PhysicalMax)) = 1e3;
    

hdr.DigitalMin = repmat(-32768,1,nchs);
hdr.DigitalMax = repmat(32767,1,nchs);




% write the edf
out_dir = [locations.main_folder,'data/eeg_clips/'];
out_name = sprintf('%s_%1.1f_to_%1.1f.edf',file_name,start_time,end_time);
edfwrite([out_dir,out_name],hdr,values,'InputSampleType',"physical");


