function select_spikes_ieeg(overwrite,do_ieeg)
    % select_spikes_gui
    % This script allows users to select epileptiform discharges in EEG data.
    % It supports an overwrite option to control processing of files.
    %
    % Usage:
    %   select_spikes_gui()           % Runs with overwrite = 0 (default)
    %   select_spikes_gui(overwrite)  % Specify overwrite as 0 or 1
    
    if nargin < 1
        overwrite = 0;  % Default to not overwriting existing results
        do_ieeg = 0;
    end

    if nargin < 2
        do_ieeg = 0;
    end
    
    clearvars -except overwrite do_ieeg

    %% File locations and set path
    locations = scalp_toolbox_locs;
    % Add script folder to path
    scripts_folder = locations.script_folder;
    addpath(genpath(scripts_folder));
    results_folder = [locations.main_folder,'results/spike_selection/'];
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end

    % ieeg stuff
    ieeg_folder = locations.ieeg_folder;
    addpath(genpath(ieeg_folder));
    pwfile = locations.ieeg_pw_file;
    login_name = locations.ieeg_login;

    % ieeg filename
    ieeg_files = {'EMU1060_Day01_1'};
    start_times = [5099.04];
    ieeg_duration = 15;
    
    % Directory containing the EDF files
    edf_dir = '/Users/erinconrad/Library/CloudStorage/Box-Box/SN12_three_threshold/SN2/threshold_2/edf_9';
    edf_files = dir(fullfile(edf_dir, '*.edf'));

    if do_ieeg
        nfiles = length(ieeg_files);
        file_names = ieeg_files;
    else
        nfiles = length(edf_files);
        file_names = {};
        for i_f = 1:nfiles
            file_names = [file_names;edf_files(i_f).name];
        end
    end
    
    % Path to the results CSV file
    results_file = fullfile(results_folder, 'spike_selections.csv');
    
    % Initialize gain and montage type
    gain = 20;
    montage_type = 'default';  % Possible values: 'default', 'b', 'c'
    
    % Check overwrite option and existing results
    if ~do_ieeg
        if overwrite == 0 && exist(results_file, 'file')
            % Load existing results
            opts = detectImportOptions(results_file, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
            existing_results = readtable(results_file, opts);
    
            processed_files = unique(existing_results.Filename);
            % Find indices of EDF files that have not been processed yet
            unprocessed_files = [];
            for i = 1:length(edf_files)
                if ~ismember(edf_files(i).name, processed_files)
                    unprocessed_files = [unprocessed_files, i];
                end
            end
            % Update edf_files to only include unprocessed files
            edf_files = edf_files(unprocessed_files);
            % Check if there are files to process
            if isempty(edf_files)
                disp('All files have been processed.');
                return;
            else
                disp(['Resuming from file: ', edf_files(1).name]);
            end
        else
            % If overwrite is 1 or results file doesn't exist, proceed as usual
            if overwrite == 1 && exist(results_file, 'file')
                delete(results_file);  % Delete existing results file
                disp('Existing results file deleted. Starting from the first file.');
            else
                disp('Starting from the first file.');
            end
            % edf_files remains as is
        end
    end
    
    % Loop through each EDF file
    for i = 1:nfiles
        if do_ieeg
            data = download_ieeg_data(file_names{i},login_name,pwfile,...
                [start_times(i) start_times(i)+ieeg_duration],1); % 1 means get lots of data
            values = data.values;
            chLabels = data.chLabels(:,1);
            fs = data.fs;

            non_intracranial = find_non_intracranial(chLabels);
            values = values(:,~non_intracranial);
            chLabels = chLabels(~non_intracranial);

        else
            % Get the full path of the current EDF file
            edf_filename = fullfile(edf_dir, edf_files(i).name);
            
            % Load and process the EDF file
            [values, chLabels, fs] = read_in_edf(edf_filename);
            data.fs = fs;
            data.values = values;
        end
        
        % Preprocess the data
        values = preprocess_eeg(values, fs);
        
        % Apply initial montage
        [processed_values, processed_labels] = apply_montage(chLabels, values, montage_type);
        
        % Initialize figure
        hFig = figure;
        set(hFig, 'Position', [10 10 1400 1000]);
        ax = axes(hFig);
        plot_handles = plot_eeg(ax, processed_values, fs, processed_labels, gain);
        hold on  % Allow adding highlights
        title(['File: ', file_names{i}], 'Interpreter', 'none');
        setappdata(hFig, 'axes_handle', ax);

        % Add a text box for displaying messages
        message_text = uicontrol('Style', 'text', 'Parent', hFig, 'Units', 'normalized', ...
            'Position', [0.1, 0.95, 0.8, 0.03], 'String', '', 'FontSize', 12, ...
            'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
        
        % Add a button for channel selection
        channel_button = uicontrol('Style', 'pushbutton', 'Parent', hFig, 'Units', 'normalized', ...
            'Position', [0.82, 0.95, 0.08, 0.04], 'String', 'Channels', 'FontSize', 12, ...
            'Callback', @channelButtonCallback);
        
        % Initialize variables for appdata
        setappdata(hFig, 'gain', gain);
        setappdata(hFig, 'montage_type', montage_type);
        setappdata(hFig, 'values', values);
        setappdata(hFig, 'chLabels', chLabels);
        setappdata(hFig, 'fs', fs);
        setappdata(hFig, 'processed_values', processed_values);
        setappdata(hFig, 'processed_labels', processed_labels);
        setappdata(hFig, 'plot_handles', plot_handles);
        setappdata(hFig, 'selected_spikes', {});  % Store selected spikes as empty cell array
        setappdata(hFig, 'highlight_patches', {});  % Store handles to highlight patches
        setappdata(hFig, 'confirmation_stage', false);  % Flag for confirmation
        setappdata(hFig, 'montage_warning_stage', false);
        setappdata(hFig, 'requested_montage', '');
        setappdata(hFig, 'message_text', message_text);  % Store message text handle
        setappdata(hFig, 'channel_visibility', []);  % Initialize channel visibility
        setappdata(hFig, 'axes_handle', ax);  % Store axes handle
        setappdata(hFig, 'edf_files', edf_files);  % Store edf_files
        setappdata(hFig, 'current_file_index', i);  % Store current file index
        
        
        % Set up figure callbacks
        set(hFig, 'WindowKeyPressFcn', @keyPressHandler);
        set(hFig, 'WindowButtonDownFcn', @mouseClickHandler);
        set(hFig, 'CloseRequestFcn', @figureCloseRequest);
        
        % Wait for user interaction
        uiwait(hFig);
        
        % Check if the figure is still valid (user may have closed it)
        if ~ishandle(hFig)
            disp('Figure closed by user. Exiting program.');
            break;  % Exit the for loop
        end
        
        % Retrieve selected spikes
        selected_spikes = getappdata(hFig, 'selected_spikes');
        fs = getappdata(hFig, 'fs');
        processed_labels = getappdata(hFig, 'processed_labels');
        
        % Close the figure
        if ishandle(hFig)
            close(hFig);
        end
        
        % Save the spike selections to the CSV file
        save_spike_selections({file_names{i}, selected_spikes}, results_file, fs, processed_labels);
        disp(['Spike selections for ', file_names{i}, ' saved.']);
    end
    
    disp('Processing completed.');
    
    %% Nested functions
    
    function values = preprocess_eeg(values, fs)
        % Preprocess EEG data: interpolate NaNs, filter, demean
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
    end
    
    function [processed_values, processed_labels] = apply_montage(chLabels, values, montage_type)
        % Apply the selected montage
        switch montage_type
            case 'default'
                % Default montage (e.g., bipolar montage)
                if do_ieeg
                    [processed_values,~,processed_labels] =...
                        bipolar_montage(values,chLabels,[],[],[],[]);
                else
                    [processed_values, processed_labels] = scalp_bipolar(chLabels, values);
                end
            case 'b'
                if do_ieeg
                    [processed_values,~,processed_labels] =...
                        bipolar_montage(values,chLabels,[],[],[],[]);
                else
                    [processed_values, processed_labels] = scalp_bipolar(chLabels, values);
                end
            case 'c'
                % Alternative montage 'c' (user-defined)
                [processed_values, processed_labels] = car_montage_2(values,chLabels);
            otherwise
                if do_ieeg
                    [processed_values,~,processed_labels] =...
                        bipolar_montage(values,chLabels,[],[],[],[]);
                else
                    [processed_values, processed_labels] = scalp_bipolar(chLabels, values);
                end
        end
    end
    
    function plot_handles = plot_eeg(ax, values, fs, labels, gain)
        % Plot the EEG data with the specified gain
        % Calculate time vector
        t = (0:size(values, 1) - 1) / fs;
        
        % Scale the values according to the gain and channel offsets
        n_channels = size(values, 2);
        offsets = (1:n_channels)' * gain * 10;  % Adjust multiplier as needed
        scaled_values = values + offsets';
        
        % Plot each channel and store handles
        plot_handles = gobjects(n_channels, 1);
        for ch = 1:n_channels
            plot_handles(ch) = plot(ax, t, scaled_values(:, ch), 'Color', 'k', 'HitTest', 'off', 'PickableParts', 'none');
            hold on
        end
        
        % Adjust axes properties
        ax.YDir = 'reverse';  % Reverse y-axis to match EEG convention
        ax.YTick = offsets;
        ax.YTickLabel = labels;
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Channels');
        ylim(ax, [0, max(offsets) + gain * 10]);
        xlim(ax, [t(1), t(end)]);
    end

    
    function keyPressHandler(hObject, event)
        % Handle key presses
        gain = getappdata(hObject, 'gain');
        montage_type = getappdata(hObject, 'montage_type');
        values = getappdata(hObject, 'values');
        chLabels = getappdata(hObject, 'chLabels');
        fs = getappdata(hObject, 'fs');
        confirmation_stage = getappdata(hObject, 'confirmation_stage');
        requested_montage = getappdata(hObject, 'requested_montage');
        montage_warning_stage = getappdata(hObject, 'montage_warning_stage');
        selected_spikes = getappdata(hObject, 'selected_spikes');
        channel_visibility = getappdata(hObject, 'channel_visibility');
        
        if strcmp(event.Key, 'uparrow')
            update_message(hObject, 'Gain decreased.');
            gain = gain - 10;
            setappdata(hObject, 'gain', gain);
            redraw_plot(hObject);
            % Reset montage warning stage
            setappdata(hObject, 'montage_warning_stage', false);
            setappdata(hObject, 'requested_montage', '');
        elseif strcmp(event.Key, 'downarrow')
            update_message(hObject, 'Gain increased.');
            gain = gain + 10;
            setappdata(hObject, 'gain', gain);
            redraw_plot(hObject);
            % Reset montage warning stage
            setappdata(hObject, 'montage_warning_stage', false);
            setappdata(hObject, 'requested_montage', '');
        elseif strcmp(event.Key, 'b') || strcmp(event.Key, 'c')
            % Montage change requested
            %requested_montage = event.Key;
            if ~isempty(selected_spikes)
                if ~montage_warning_stage || ~strcmp(requested_montage, event.Key)
                    update_message(hObject, ['Changing the montage will de-select all currently selected spikes. ', ...
                        'Press "', event.Key, '" again to confirm montage change.']);
                    setappdata(hObject, 'montage_warning_stage', true);
                    setappdata(hObject, 'requested_montage', requested_montage);
                else
                    % User confirms montage change
                    update_message(hObject, ['Montage changed to ', event.Key, '. Spikes de-selected.']);
                    montage_type = event.Key;
                    setappdata(hObject, 'montage_type', montage_type);
                    setappdata(hObject, 'selected_spikes', {});  % Clear selected spikes
                    setappdata(hObject, 'highlight_patches', {});  % Clear highlights
                    setappdata(hObject, 'channel_visibility', []);  % Reset channel visibility
                    setappdata(hObject, 'montage_warning_stage', false);
                    setappdata(hObject, 'requested_montage', '');
                    redraw_plot(hObject);
                end
            else
                % No spikes selected, proceed with montage change
                update_message(hObject, ['Montage changed to ', event.Key]);
                montage_type = event.Key;
                setappdata(hObject, 'montage_type', montage_type);
                setappdata(hObject, 'channel_visibility', []);  % Reset channel visibility

                redraw_plot(hObject);
            end
        elseif strcmp(event.Key, 'return')
            if ~confirmation_stage
                update_message(hObject, 'Press Enter again to confirm spike selections and proceed.');
                setappdata(hObject, 'confirmation_stage', true);
            else
                update_message(hObject, 'Spike selections confirmed.');
                uiresume(hObject);
            end
            % Reset montage warning stage
            setappdata(hObject, 'montage_warning_stage', false);
            setappdata(hObject, 'requested_montage', '');
        else
            update_message(hObject, ['Key pressed: ', event.Key]);
            % Reset montage warning stage
            setappdata(hObject, 'montage_warning_stage', false);
            setappdata(hObject, 'requested_montage', '');
        end
    end

    
   
    function mouseClickHandler(hObject, ~)
    % Handle mouse clicks
    ax = hObject.CurrentAxes;
    fs = getappdata(hObject, 'fs');
    processed_values = getappdata(hObject, 'processed_values');
    processed_labels = getappdata(hObject, 'processed_labels');
    selected_spikes = getappdata(hObject, 'selected_spikes');
    highlight_patches = getappdata(hObject, 'highlight_patches');
    gain = getappdata(hObject, 'gain');
    channel_visibility = getappdata(hObject, 'channel_visibility');

    % Get handles to plot and labels for highlighting
    plot_handles = getappdata(hObject, 'plot_handles');  % Handles for plotted EEG channels
    ytick_labels = ax.YTickLabel;  % Get the current y-tick labels (channel labels)

    % Handle temporary highlighting
    highlighted_trace = getappdata(hObject, 'highlighted_trace');
    highlighted_label_idx = getappdata(hObject, 'highlighted_label_idx');

    % Ensure selected_spikes is a cell array
    if isempty(selected_spikes)
        selected_spikes = {};
        highlight_patches = {};
    end

    % Filter channels based on visibility
    if isempty(channel_visibility)
        channel_visibility = true(length(processed_labels), 1);
    end
    visible_channels = find(channel_visibility);
    filtered_labels = processed_labels(visible_channels);
    filtered_values = processed_values(:, visible_channels);

    % Get click coordinates
    clickPoint = get(ax, 'CurrentPoint');
    x_click = clickPoint(1, 1);  % Time
    y_click = clickPoint(1, 2);  % Channel offset

    % Find the closest channel in filtered list
    n_channels = length(filtered_labels);
    offsets = (1:n_channels)' * gain * 10;  % Same offset calculation
    [~, visible_ch_idx] = min(abs(offsets - y_click));
    ch_idx = visible_channels(visible_ch_idx);  % Map back to original channel index
    channel = processed_labels{ch_idx};

    % Find the closest time sample
    time = x_click;
    sample_idx = round(time * fs);
    sample_idx = max(min(sample_idx, size(processed_values, 1)), 1);  % Ensure within bounds

    % Define sample tolerance for de-selection (±0.5 seconds)
    sample_tolerance = round(0.5 * fs);

    % Check if this spike is already selected (within tolerance)
    existing_idx = find(cellfun(@(x) x(1) == ch_idx && abs(x(2) - sample_idx) <= sample_tolerance, selected_spikes), 1);

    if isempty(existing_idx)
        % Add spike to selections
        selected_spikes{end+1} = [ch_idx, sample_idx];

        % Highlight the selected area if channel is visible
        if channel_visibility(ch_idx)
            % Define time window around the click (±0.25 seconds for highlight)
            window_samples = round(0.25 * fs);  % Half a second total for highlight
            start_idx = max(sample_idx - window_samples, 1);
            end_idx = min(sample_idx + window_samples, size(processed_values, 1));
            t = (start_idx:end_idx) / fs;

            % Get channel offset in filtered channels
            n_channels = length(filtered_labels);
            offsets = (1:n_channels)' * gain * 10;
            y_offset = offsets(visible_ch_idx);
            y_data = [y_offset - gain*5, y_offset + gain*5];
            hPatch = patch(ax, [t, fliplr(t)], [repmat(y_data(1), 1, length(t)), repmat(y_data(2), 1, length(t))], ...
                'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            highlight_patches{end+1} = hPatch;
        else
            % Channel is not visible, so we don't draw the highlight
            highlight_patches{end+1} = [];
        end
        % Update message
        update_message(hObject, ['Spike selected at ', processed_labels{ch_idx}, ', Time: ', num2str(time, '%.2f'), ' s']);
    else
        % Remove spike from selections
        selected_spikes(existing_idx) = [];
        % Remove highlight if it exists
        if ~isempty(highlight_patches{existing_idx})
            delete(highlight_patches{existing_idx});
        end
        highlight_patches(existing_idx) = [];
        % Update message
        update_message(hObject, ['Spike de-selected at ', processed_labels{ch_idx}, ', Time: ', num2str(time, '%.2f'), ' s']);
    end

    % --- Temporary highlight of most recent selected spike ---

    % Restore previous highlight (if any)
    if ~isempty(highlighted_trace)
        set(highlighted_trace, 'Color', 'k', 'LineWidth', 1);  % Reset previous trace to black and normal thickness
    end
    if ~isempty(highlighted_label_idx)
        set(ax, 'FontWeight', 'normal');  % Reset previous label to normal
    end

    % Highlight the current selected spike trace and label
    set(plot_handles(visible_ch_idx), 'Color', 'b', 'LineWidth', 2);  % Highlight in blue
    set(ax, 'FontWeight', 'bold');  % Make current label bold

    % Store the new highlighted trace and label index
    setappdata(hObject, 'highlighted_trace', plot_handles(visible_ch_idx));
    setappdata(hObject, 'highlighted_label_idx', visible_ch_idx);

    % Update appdata
    setappdata(hObject, 'selected_spikes', selected_spikes);
    setappdata(hObject, 'highlight_patches', highlight_patches);
end



  
    
    function redraw_plot(hObject)
        % Redraw the EEG plot with current settings
        gain = getappdata(hObject, 'gain');
        montage_type = getappdata(hObject, 'montage_type');
        values = getappdata(hObject, 'values');
        chLabels = getappdata(hObject, 'chLabels');
        fs = getappdata(hObject, 'fs');
        selected_spikes = getappdata(hObject, 'selected_spikes');
        highlight_patches = getappdata(hObject, 'highlight_patches');
        ax = getappdata(hObject, 'axes_handle');  % Retrieve axes handle
        channel_visibility = getappdata(hObject, 'channel_visibility');

        
        % Apply montage
        [processed_values, processed_labels] = apply_montage(chLabels, values, montage_type);
        setappdata(hObject, 'processed_values', processed_values);
        setappdata(hObject, 'processed_labels', processed_labels);

         % If channel_visibility is empty, default to showing all channels
        if isempty(channel_visibility)
            channel_visibility = true(length(processed_labels), 1);
        end
        % Update channel_visibility in appdata
        setappdata(hObject, 'channel_visibility', channel_visibility);
    
        % Filter channels based on visibility
        visible_channels = find(channel_visibility);
        filtered_values = processed_values(:, visible_channels);
        filtered_labels = processed_labels(visible_channels);
        
        % Clear axes and re-plot
        cla(ax); 
        plot_eeg(ax, filtered_values, fs, filtered_labels, gain);
        hold(ax,'on');
        
        % Redraw highlights
        % Initialize a new cell array for highlight patches
        new_highlight_patches = cell(size(selected_spikes));
        for idx = 1:length(selected_spikes)
            spike_id = selected_spikes{idx};
            ch_idx = spike_id(1);
            sample_idx = spike_id(2);

            % Check if this channel is visible
            if channel_visibility(ch_idx)
                % Find the index of this channel in the filtered list
                visible_ch_idx = find(visible_channels == ch_idx);

                % Define time window around the spike
                window_samples = round(0.25 * fs);  % Half a second total
                start_idx = max(sample_idx - window_samples, 1);
                end_idx = min(sample_idx + window_samples, size(processed_values, 1));
                t = (start_idx:end_idx) / fs;

                % Get channel offset
                n_channels = length(filtered_labels);
                offsets = (1:n_channels)' * gain * 10;
                y_offset = offsets(visible_ch_idx);
                y_data = [y_offset - gain*5, y_offset + gain*5];
                hPatch = patch(ax, [t, fliplr(t)], [repmat(y_data(1), 1, numel(t)), repmat(y_data(2), 1, numel(t))], ...
                    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                new_highlight_patches{idx} = hPatch;
            else
                % Channel is not visible, skip highlighting
                new_highlight_patches{idx} = [];
            end
        end
        % Update highlight_patches in appdata
        setappdata(hObject, 'highlight_patches', new_highlight_patches);

        
        title(['File: ', file_names{i}], 'Interpreter', 'none');
    end
    
    function update_message(hObject, message)
        % Update the message displayed in the GUI
        message_text = getappdata(hObject, 'message_text');
        set(message_text, 'String', message);
    end

    function channelButtonCallback(hObject, ~)
    % Callback function for the Channels button
    hFig = ancestor(hObject, 'figure');
    processed_labels = getappdata(hFig, 'processed_labels');
    channel_visibility = getappdata(hFig, 'channel_visibility');

    % Prepare data for the table
    num_channels = length(processed_labels);
    if isempty(channel_visibility)
        channel_visibility = true(num_channels, 1);
    end

    % Create the figure for channel selection
    dlg_title = 'Select Channels to Display';
    selectionFig = figure('Name', dlg_title, 'NumberTitle', 'off', 'MenuBar', 'none', ...
        'ToolBar', 'none', 'Resize', 'off', 'WindowStyle', 'modal');
    set(selectionFig, 'Position', [100, 100, 300, 400]);

    % Prepare data for the table
    data = [num2cell(channel_visibility), processed_labels];
    % Column names
    colnames = {'Show', 'Channel'};
    % Column format
    columnformat = {'logical', 'char'};
    % Column editable
    columneditable = [true, false];

    % Create the uitable
    t = uitable('Parent', selectionFig, 'Units', 'normalized', 'Position', [0.05, 0.15, 0.9, 0.8], ...
        'Data', data, 'ColumnName', colnames, 'ColumnFormat', columnformat, ...
        'ColumnEditable', columneditable, 'RowName', [], 'CellSelectionCallback', @cellSelectionCallback);

    % Store selected rows in appdata
    setappdata(selectionFig, 'selectedRows', []);

    % Add Select, Deselect, OK, and Cancel buttons
    select_button = uicontrol('Style', 'pushbutton', 'Parent', selectionFig, 'Units', 'normalized', ...
        'Position', [0.05, 0.02, 0.2, 0.1], 'String', 'Select', 'FontSize', 12, ...
        'Callback', @selectButtonCallback);

    deselect_button = uicontrol('Style', 'pushbutton', 'Parent', selectionFig, 'Units', 'normalized', ...
        'Position', [0.28, 0.02, 0.2, 0.1], 'String', 'Deselect', 'FontSize', 12, ...
        'Callback', @deselectButtonCallback);

    ok_button = uicontrol('Style', 'pushbutton', 'Parent', selectionFig, 'Units', 'normalized', ...
        'Position', [0.51, 0.02, 0.2, 0.1], 'String', 'OK', 'FontSize', 12, ...
        'Callback', @okButtonCallback);

    cancel_button = uicontrol('Style', 'pushbutton', 'Parent', selectionFig, 'Units', 'normalized', ...
        'Position', [0.74, 0.02, 0.2, 0.1], 'String', 'Cancel', 'FontSize', 12, ...
        'Callback', @cancelButtonCallback);

    % Callback functions
    function cellSelectionCallback(~, event)
        % Keep track of selected rows, ensuring no duplicates
        selectedRows = unique(event.Indices(:,1));
        setappdata(selectionFig, 'selectedRows', selectedRows);
    end

    function selectButtonCallback(~, ~)
        % Toggle 'Select' for the selected rows
        selectedRows = getappdata(selectionFig, 'selectedRows');
        if isempty(selectedRows)
            return;
        end
        data = get(t, 'Data');
        for idx = 1:length(selectedRows)
            data{selectedRows(idx), 1} = true;  % Set checkboxes to true (select)
        end
        set(t, 'Data', data);  % Update the table
    end

    function deselectButtonCallback(~, ~)
        % Toggle 'Deselect' for the selected rows
        selectedRows = getappdata(selectionFig, 'selectedRows');
        if isempty(selectedRows)
            return;
        end
        data = get(t, 'Data');
        for idx = 1:length(selectedRows)
            data{selectedRows(idx), 1} = false;  % Set checkboxes to false (deselect)
        end
        set(t, 'Data', data);  % Update the table
    end

    function okButtonCallback(~, ~)
        % Get the data from the table
        table_data = get(t, 'Data');
        new_visibility = cell2mat(table_data(:, 1));
        % Update the channel visibility in appdata
        setappdata(hFig, 'channel_visibility', new_visibility);
        % Close the selection figure
        close(selectionFig);
        % Redraw the plot with updated visibility
        redraw_plot(hFig);
    end

    function cancelButtonCallback(~, ~)
        % Close the selection figure without changing anything
        close(selectionFig);
    end
end


    function figureCloseRequest(hObject, ~)
        disp('Figure close requested by user.');
        % Resume execution in case uiwait is blocking
        uiresume(hObject);
        % Delete the figure
        delete(hObject);
    end
    
    function save_spike_selections(selection, filename, fs, processed_labels)
        % Save the spike selections to a CSV file
        % The CSV file will have columns: Filename,Channel,Time_s
        % 'selection' is a cell array: {filename_i, selected_spikes}
        filename_i = selection{1};
        spikes = selection{2};
    
        % Check if the file exists to determine whether to write header
        write_header = ~exist(filename, 'file');
    
        fid = fopen(filename, 'a');  % Open in append mode
        if write_header
            % Write header with valid variable names
            fprintf(fid, 'Filename,Channel,Time_s\n');
        end
        if isempty(spikes)
            % Write an entry indicating no spikes selected
            fprintf(fid, '%s,None,NaN\n', filename_i);
        else
            for s = 1:length(spikes)
                ch_idx = spikes{s}(1);
                sample_idx = spikes{s}(2);
                time_s = sample_idx / fs;
                fprintf(fid, '%s,%s,%.3f\n', filename_i, processed_labels{ch_idx}, time_s);
            end
        end
        fclose(fid);
    end


    
  

end
