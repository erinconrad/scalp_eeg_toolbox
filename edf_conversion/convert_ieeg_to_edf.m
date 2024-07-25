%{
file_name = 'EMU2099_Day01_1';
start_time = 12160.30;
end_time = 12160.30+15;
%}
file_name = 'EMU1732_Day01_1';
start_time = 5266.42;
end_time = 5285.32;

%% File locs and set path
locations = scalp_toolbox_locs;
% add script folder to path
scripts_folder = locations.script_folder;
edf_files = [locations.main_folder,'results/edf_files/'];
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

%% Make edf header
hdr = make_header_edf(values,chLabels,fs);

%% Convert to edf and save
out_file = sprintf('%s_%1.1fs_%1.1fs.edf',file_name,start_time,end_time);
edfw = edfwrite([edf_files,out_file],hdr,values,'InputSampleType',"physical");