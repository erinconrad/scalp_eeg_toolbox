% validate spike detections

clear

%% overwrite settings
overwrite = 2;

%% Name output file and edf directory
edf_dir = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/example_edfs/';
subfolder_name = 'edf_Nov1924';
out_file = 'persyst_50_test.csv';

%% File locations and set path
locations = scalp_toolbox_locs;
% Add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
results_folder = [locations.main_folder,'results/spike_net_validation/'];

%% Define variable names in table
% Define the variable names
variableNames = {'Foldername', 'Filename', 'Response'};  % Replace with your desired variable names

if overwrite == 0
    % Load the out file if it exists
    if exist([results_folder,out_file],"file") ~= 0
        responses = readtable([results_folder,out_file],...
            'ReadVariableNames', true,'HeaderLines', 0);
    else
        responses = table('Size', [0, 3], 'VariableTypes', {'string', 'string', 'string'}, ...
                  'VariableNames', variableNames);
    end
else
    responses = table('Size', [0, 3], 'VariableTypes', {'string', 'string', 'string'}, ...
                  'VariableNames', variableNames);
end



%% Get subfolders
% Get the directory contents
contents = dir(edf_dir);

% Filter out the subfolders (ignoring '.' and '..')
subfolders = contents([contents.isdir]);  % Get only directories
subfolderNames = {subfolders.name};  % Extract the names
subfolderNames = subfolderNames(~ismember(subfolderNames, {'.', '..'}));  % Remove '.' and '..'

% First, loop over subfolders
for s = 1:length(subfolderNames)

    if ~isempty(subfolder_name)
        if ~strcmp(subfolderNames{s},subfolder_name)
            continue
        end
    end

    sub_dir = [edf_dir,subfolderNames{s},'/'];

    % Get the list of edf files
    edf_files = dir(fullfile(sub_dir, '*.edf'));


    % Loop through each EDF file
    for i = 1:length(edf_files)
        % Get the full path of the current EDF file
        edf_filename = fullfile(sub_dir, edf_files(i).name);

        % Look and see if we've done it already (if not overwriting)
        if overwrite == 0 || overwrite == 2
            completed_folders = responses.Foldername;
            completed_files = responses.Filename;

            if any(strcmp(completed_folders,subfolderNames{s}) & ...
                    strcmp(completed_files,edf_files(i).name))
                fprintf('\nAlready did %s and %s, skipping\n',...
                    subfolderNames{s},edf_files(i).name);
                continue;
            end
        end

        % Initialize gain
        gain = 20;
        
        % Load and process the EDF file
        [values, chLabels, fs] = read_in_edf(edf_filename);
        data.fs = fs;
        data.values = values;
        
        % Interpolate over NaNs
        for j = 1:size(values, 2)
            nan_indices = isnan(values(:, j));
            if any(nan_indices)
                values(nan_indices, j) = nanmean(values(:, j));
            end
        end

        % Demean the channels
        values = values - nanmean(values, 1);
        
        % High-pass filter
        values = highpass(values, 1, fs);
        values = bandstop(values, [58 62], fs);
        
        
        
        [bipolar_values, bipolar_labels] = scalp_bipolar(chLabels, values);
        data.bipolar_values = bipolar_values;
        data.bipolar_labels = bipolar_labels;

        % also get car
        [car_values,car_labels] = car_montage_2(values,chLabels);
        
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

        % Create a figure for displaying the EEG data
        figure
        hFig = gcf;
        set(hFig, 'Position', [10 10 1400 1000]);
        tiledlayout(1, 1, 'Padding', 'compact');
        nexttile
        plot_scalp_eeg_gain(plot_values, fs, plot_labels, gain);
        fprintf('\nDisplaying %s\n',(edf_files(i).name))
        
        % Initialize variables and store them in appdata
        user_input = '';
        arrow_key_input = '';
        setappdata(hFig, 'user_input', user_input);
        setappdata(hFig, 'arrow_key_input', arrow_key_input);
        setappdata(hFig, 'gain', gain);
        setappdata(hFig, 'bipolar_values', bipolar_values);
        setappdata(hFig, 'bipolar_labels', bipolar_labels);
        setappdata(hFig, 'car_values', car_values);
        setappdata(hFig, 'car_labels', car_labels);
        setappdata(hFig, 'fs', fs);
        setappdata(hFig, 'montage_type', montage_type);
        
        % Set up the figure to detect key presses using WindowKeyPressFcn
        set(hFig, 'WindowKeyPressFcn', @keyPressHandler);
        
        % Wait until the user has confirmed their response
        uiwait(hFig);  % Wait until uiresume is called
        
        % Retrieve the stored arrow_key_input and user_input from appdata
        arrow_key_input = getappdata(hFig, 'arrow_key_input');
        user_input = getappdata(hFig, 'user_input');
        gain = getappdata(hFig, 'gain');
        montage_type = getappdata(hFig, 'montage_type');
        
        % Now you can safely close the figure
        if ishandle(hFig)
            close(hFig);
        end
        
        % Store the filename and response in the responses array
        responses(end+1,:) = {subfolderNames{s}, edf_files(i).name, user_input};
        
        if overwrite ~=2
            % Save the table as a CSV file (save as you go)
            writetable(responses, [results_folder, out_file]);
        end

    end
end



% Nested key press handler function
function keyPressHandler(hObject, event)
    % Retrieve data stored in appdata
    arrow_key_input = getappdata(hObject, 'arrow_key_input');
    user_input = getappdata(hObject, 'user_input');
    gain = getappdata(hObject, 'gain');
    bipolar_values = getappdata(hObject, 'bipolar_values');
    bipolar_labels = getappdata(hObject, 'bipolar_labels');
    car_values = getappdata(hObject, 'car_values');
    car_labels = getappdata(hObject, 'car_labels');
    fs = getappdata(hObject, 'fs');
    montage_type = getappdata(hObject, 'montage_type');

    % Handle key presses
    if strcmp(event.Key, 'uparrow')
        disp('Up arrow pressed');
        arrow_key_input = 'up';
        gain = gain - 10;
        setappdata(hObject, 'gain', gain);
        hold off
        % Re-plot with new montage
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        plot_scalp_eeg_gain(plot_values, fs, plot_labels, gain);
        
    elseif strcmp(event.Key, 'downarrow')
        disp('Down arrow pressed');
        arrow_key_input = 'down';
        gain = gain + 10;
        setappdata(hObject, 'gain', gain);
        hold off
        % Re-plot with new montage
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        plot_scalp_eeg_gain(plot_values, fs, plot_labels, gain);
    elseif strcmp(event.Key, 'b') || strcmp(event.Key, 'c')
            disp(['Montage change requested: ', event.Key]);
            if strcmp(event.Key,'b')
                montage_type = 'bipolar';
            elseif strcmp(event.Key,'c')
                montage_type = 'car';
            end
            setappdata(hObject, 'montage_type', montage_type);
            hold off
            % Re-plot with new montage
            switch montage_type
                case 'bipolar'
                    plot_values = bipolar_values;
                    plot_labels = bipolar_labels;
                case 'car'
                    plot_values = car_values;
                    plot_labels = car_labels;
            end
            plot_scalp_eeg_gain(plot_values, fs, plot_labels, gain);
    elseif strcmp(event.Key, 'return')
        % User pressed Enter to confirm input
        if strcmp(user_input, 'y') || strcmp(user_input, 'n') || strcmp(user_input, 'f')
            disp(['User input confirmed: ', user_input]);
            % Store the user input
            setappdata(hObject, 'user_input', user_input);
            % Resume execution (but do NOT close the figure here)
            uiresume(hObject);
        else
            disp('Please press "y" or "n" before pressing Enter to confirm.');
        end
    elseif strcmp(event.Key, 'y') || strcmp(event.Key, 'n') || strcmp(event.Key, 'f')
        % Store the user's input character but wait for Enter to confirm
        user_input = event.Key;
        disp(['You pressed "', user_input, '". Press Enter to confirm.']);
        % Store the tentative user input
        setappdata(hObject, 'user_input', user_input);
    else
        disp(['Other key pressed: ', event.Key]);
    end

    % Store the updated data back into appdata
    setappdata(hObject, 'arrow_key_input', arrow_key_input);
end
