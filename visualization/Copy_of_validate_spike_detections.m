function validate_spike_detections()
    % validate spike detections

    clear

    %% File locations and set path
    locations = scalp_toolbox_locs;
    % Add script folder to path
    scripts_folder = locations.script_folder;
    addpath(genpath(scripts_folder));
    results_folder = [locations.main_folder,'results/spike_net_validation/sn1/thresh_0.8/'];

    edf_dir = '/Users/erinconrad/Library/CloudStorage/Box-Box/SN12_three_threshold/SN1/edf_threshold_8/';
    edf_files = dir(fullfile(edf_dir, '*.edf'));

    % Initialize an empty cell array to store filenames and user responses
    responses = {}; 

    % Initialize gain
    gain = 20;

    % Loop through each EDF file
    for i = 1:length(edf_files)
        % Get the full path of the current EDF file
        edf_filename = fullfile(edf_dir, edf_files(i).name);

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

        % High-pass filter
        values = highpass(values, 0.5, fs);
        values = bandstop(values, [58 62], fs);

        % Demean the channels
        values = values - nanmean(values, 1);

        [bipolar_values, bipolar_labels] = scalp_bipolar(chLabels, values);
        data.bipolar_values = bipolar_values;
        data.bipolar_labels = bipolar_labels;

        % Create a figure for displaying the EEG data
        figure
        hFig = gcf;
        set(hFig, 'Position', [10 10 1400 1000]);
        tiledlayout(1, 1, 'Padding', 'compact');
        nexttile
        plot_scalp_eeg_gain(bipolar_values, fs, bipolar_labels, gain);

        % Initialize variables and store them in appdata
        user_input = '';
        arrow_key_input = '';
        setappdata(hFig, 'user_input', user_input);
        setappdata(hFig, 'arrow_key_input', arrow_key_input);
        setappdata(hFig, 'gain', gain);
        setappdata(hFig, 'bipolar_values', bipolar_values);
        setappdata(hFig, 'bipolar_labels', bipolar_labels);
        setappdata(hFig, 'fs', fs);

        % Set up the figure to detect key presses using WindowKeyPressFcn
        set(hFig, 'WindowKeyPressFcn', @keyPressHandler);

        % Wait until the user has confirmed their response
        uiwait(hFig);  % Wait until uiresume is called

        % Retrieve the stored arrow_key_input and user_input from appdata
        arrow_key_input = getappdata(hFig, 'arrow_key_input');
        user_input = getappdata(hFig, 'user_input');
        gain = getappdata(hFig, 'gain');

        % Now you can safely close the figure
        if ishandle(hFig)
            close(hFig);
        end

        % Store the filename and response in the responses array
        responses = [responses; {edf_files(i).name, user_input}];

        % Optionally, you can save the gain if needed
        % gains{i} = gain;
    end

    % Convert the responses cell array into a table
    response_table = cell2table(responses, 'VariableNames', {'Filename', 'Response'});

    % Save the table as a CSV file
    writetable(response_table, [results_folder, 'edf_responses.csv']);
    
    disp('Responses saved to edf_responses.csv');

    % Nested key press handler function
    function keyPressHandler(hObject, event)
        % Retrieve data stored in appdata
        arrow_key_input = getappdata(hObject, 'arrow_key_input');
        user_input = getappdata(hObject, 'user_input');
        gain = getappdata(hObject, 'gain');
        bipolar_values = getappdata(hObject, 'bipolar_values');
        bipolar_labels = getappdata(hObject, 'bipolar_labels');
        fs = getappdata(hObject, 'fs');

        % Handle key presses
        if strcmp(event.Key, 'uparrow')
            disp('Up arrow pressed');
            arrow_key_input = 'up';
            gain = gain - 10;
            setappdata(hObject, 'gain', gain);
            hold off
            % Re-plot with updated gain
            plot_scalp_eeg_gain(bipolar_values, fs, bipolar_labels, gain);
        elseif strcmp(event.Key, 'downarrow')
            disp('Down arrow pressed');
            arrow_key_input = 'down';
            gain = gain + 10;
            setappdata(hObject, 'gain', gain);
            hold off
            % Re-plot with updated gain
            plot_scalp_eeg_gain(bipolar_values, fs, bipolar_labels, gain);
        elseif strcmp(event.Key, 'return')
            % User pressed Enter to confirm input
            if strcmp(user_input, 'y') || strcmp(user_input, 'n')
                disp(['User input confirmed: ', user_input]);
                % Store the user input
                setappdata(hObject, 'user_input', user_input);
                % Resume execution (but do NOT close the figure here)
                uiresume(hObject);
            else
                disp('Please press "y" or "n" before pressing Enter to confirm.');
            end
        elseif strcmp(event.Key, 'y') || strcmp(event.Key, 'n')
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
end
