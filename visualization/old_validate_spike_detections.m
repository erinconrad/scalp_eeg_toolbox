%old_validate_spike_detections

% validate spike detections

clear

%% File locs and set path
locations = scalp_toolbox_locs;
% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
results_folder = [locations.main_folder,'results/spike_net_validation/sn2/thresh_0.8/'];

edf_dir = '/Users/erinconrad/Library/CloudStorage/Box-Box/SN12_three_threshold/SN2/edf_threshold_8/';
edf_files = dir(fullfile(edf_dir, '*.edf'));

% Initialize an empty cell array to store filenames and user responses
responses = {}; 

% Loop through each EDF file
for i = 1:length(edf_files)
    % Get the full path of the current EDF file
    edf_filename = fullfile(edf_dir, edf_files(i).name);
    
    % load the edf file
    [values,chLabels,fs] = read_in_edf(edf_filename);
    data.fs = fs;
    data.values = values;

    % Interpolate over nans
    for j = 1:size(values,2)
        values(isnan(values(:,j)),j) = nanmean(values(:,j));
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
    plot_scalp_eeg(bipolar_values,fs,bipolar_labels);

    % Prompt the user for a response (y or n)
    valid_input = false;
    while ~valid_input
        user_input = input('Enter "y" if spike or "n" if not: ', 's');
        if strcmp(user_input, 'y') || strcmp(user_input, 'n')
            valid_input = true;
        else
            disp('Invalid input. Please enter "y" or "n".');
        end
    end
    
    % Store the filename and response in the responses array
    responses = [responses; {edf_files(i).name, user_input}];
end

% Convert the responses cell array into a table
response_table = cell2table(responses, 'VariableNames', {'Filename', 'Response'});

% Save the table as a CSV file
writetable(response_table, [results_folder,'edf_responses.csv']);