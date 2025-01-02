function createApp()
    % Create the figure
    fig = uifigure('Position', [100 100 500 250], 'Name', 'Data Processor');

    % Input Path Components
    uilabel(fig, 'Position', [20 180 100 22], 'Text', 'Input File:');
    inputField = uieditfield(fig, 'text', 'Position', [130 180 250 22], 'Editable', 'off');
    inputButton = uibutton(fig, 'push', 'Position', [400 180 60 22], 'Text', 'Browse', ...
        'ButtonPushedFcn', @(btn,event) browseInput());

    % Output Path Components
    uilabel(fig, 'Position', [20 130 100 22], 'Text', 'Output File:');
    outputField = uieditfield(fig, 'text', 'Position', [130 130 250 22], 'Editable', 'off');
    outputButton = uibutton(fig, 'push', 'Position', [400 130 60 22], 'Text', 'Browse', ...
        'ButtonPushedFcn', @(btn,event) browseOutput());

    % Run Button
    runButton = uibutton(fig, 'push', 'Position', [200 50 100 30], 'Text', 'Run', ...
        'ButtonPushedFcn', @(btn,event) runProcessing());

    % Status Label
    statusLabel = uilabel(fig, 'Position', [20 20 460 22], 'Text', 'Status: Waiting for input...', 'HorizontalAlignment', 'left');

    % Nested functions for button callbacks
    function browseInput()
        [file, path] = uigetfile({'*.*', 'All Files'}, 'Select Input File');
        if isequal(file,0)
            return;
        else
            inputField.Value = fullfile(path, file);
        end
    end

    function browseOutput()
        [file, path] = uiputfile({'*.*', 'All Files'}, 'Select Output File');
        if isequal(file,0)
            return;
        else
            outputField.Value = fullfile(path, file);
        end
    end

    
    function processedData = processData(data)
        % Placeholder for your data processing logic
        processedData = data; % Modify as needed
    end
end

function compilable_spike_validator(varargin)

%% Paths - edit this!!!
% Default paths
inputPath = '';

% Parse input arguments
p = inputParser;
addParameter(p, 'input', '', @ischar);
parse(p, varargin{:});
inputPath = p.Results.input;
% Validate paths
if isempty(inputPath)
    error('Input path must be specified.');
end
edf_dir = inputPath;

%edf_dir = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/example_edfs/'; % Update this path
out_file = 'selections.csv'; % Output CSV file name

%% Overwrite settings - edit this as needed
% overwrite = 0 (default): skip previously completed files
% overwrite = 1: re-write selections of all files
% overwrite = 2: show all files, but do not save any selections
overwrite = 0;

%{ ----- NO NEED TO EDIT BELOW THIS ------------ %}

%% Define variable names in table
variableNames = {'Filename', 'Response'};  % Replace with your desired variable names

%% Load output file if it exists and initialize table
if overwrite == 0
    % Load the out file if it exists
    if exist(fullfile(edf_dir, out_file), "file") ~= 0
        responses = readtable(fullfile(edf_dir, out_file),...
            'ReadVariableNames', true, 'HeaderLines', 0);
    else
        responses = table('Size', [0, 2], 'VariableTypes', {'string', 'string'}, ...
                  'VariableNames', variableNames);
    end
else
    responses = table('Size', [0, 2], 'VariableTypes', {'string', 'string'}, ...
                  'VariableNames', variableNames);
end

% Get the list of EDF files
edf_files = dir(fullfile(edf_dir, '*.edf'));

% Loop through each EDF file
for i = 1:length(edf_files)
    
    %% Figure out the file path and whether to run it
    % Get the full path of the current EDF file
    edf_filename = fullfile(edf_dir, edf_files(i).name);
    
    % Look and see if we've done it already (if not overwriting)
    if overwrite == 0 || overwrite == 2
        completed_files = responses.Filename;
        
        if any(strcmp(completed_files, edf_files(i).name))
            fprintf('\nAlready did %s, skipping\n', edf_files(i).name);
            continue;
        end
    end
    
    %% Load and process the EDF file
    [values, chLabels, fs] = read_in_edf_clean(edf_filename); % Ensure this function is available
    
    % Initialize gain
    gain = 50;
    
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
    
    % Get bipolar montage
    [bipolar_values, bipolar_labels] = scalp_bipolar_clean(values, chLabels);
    
    % Also get CAR montage
    [car_values, car_labels] = car_montage_clean(values, chLabels);
    
    %% Plotting
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
    figure;
    hFig = gcf;
    set(hFig, 'Position', [10 10 1400 1000]);
    tiledlayout(1, 1, 'Padding', 'compact');
    nexttile;
    
    % Initialize user_input
    user_input = '';  % Initialize as empty
    
    % Plot EEG data with prompt and selection texts
    [prompt_text, selection_text] = plot_eeg_with_text(plot_values, fs, plot_labels, gain, user_input);
    
    fprintf('\nDisplaying %s\n', edf_files(i).name);
    
    %% Store Handles and Data in appdata
    setappdata(hFig, 'prompt_text', prompt_text);
    setappdata(hFig, 'selection_text', selection_text);
    setappdata(hFig, 'gain', gain);
    setappdata(hFig, 'bipolar_values', bipolar_values);
    setappdata(hFig, 'bipolar_labels', bipolar_labels);
    setappdata(hFig, 'car_values', car_values);
    setappdata(hFig, 'car_labels', car_labels);
    setappdata(hFig, 'fs', fs);
    setappdata(hFig, 'montage_type', montage_type);
    
    %% Handle User Input
    % Initialize variables and store them in appdata
    arrow_key_input = '';
    setappdata(hFig, 'user_input', user_input);
    setappdata(hFig, 'arrow_key_input', arrow_key_input);
    
    % Set up the figure to detect key presses using WindowKeyPressFcn
    set(hFig, 'WindowKeyPressFcn', @keyPressHandler);
    
    % Wait until the user has confirmed their response
    uiwait(hFig);  % Wait until uiresume is called
    
    %% Retrieve User Input and Other Data
    % Retrieve the stored arrow_key_input and user_input from appdata
    arrow_key_input = getappdata(hFig, 'arrow_key_input');
    user_input = getappdata(hFig, 'user_input');
    gain = getappdata(hFig, 'gain');
    montage_type = getappdata(hFig, 'montage_type');
    
    %% Handle End of File and Write Responses
    % Now you can safely close the figure
    if ishandle(hFig)
        close(hFig);
    end
    
    % Store the filename and response in the responses array
    responses(end+1, :) = {edf_files(i).name, user_input};
    
    if overwrite ~= 2
        % Save the table as a CSV file (save as you go)
        writetable(responses, fullfile(edf_dir, out_file));
    end
    
end
end

%% Function Definitions
%function clean_temple


function [out, chLabels, fs] = read_in_edf_clean(filepath)

data = edfread(filepath);
info = edfinfo(filepath);
fs = info.NumSamples(1);

chs = info.SignalLabels;
nchs = length(chs);
chs = replace(chs, "-", "_");
chs = replace(chs, " ", "");

nsamples = size(data,1) * fs;

out = nan(nsamples,nchs);

for i = 1:nchs
    if iscell(data.(chs(i)))
        out(:,i) = cell2mat(data.(chs(i)));
    else
        continue
    end

end

chLabels = cellstr(chs);

end

function [prompt_text, selection_text] = plot_eeg_with_text(plot_values, fs, plot_labels, gain, user_input)
    % This function plots the EEG data and adds prompt and selection texts
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
    
    % Add prompt text as super title
    prompt_str = 'Select spike: l (left), r (right), m (midline), n (no spike)';
    sgtitle(prompt_str, 'FontSize', 16, 'Color', 'blue');
    
    % Add selection display text at a fixed position (e.g., top-left corner)
    if isempty(user_input)
        selection_str = 'Selected: None';
    else
        selection_str = ['Selected: ', user_input];
    end
    % Use normalized units to position the text relative to the axes
    selection_text = text(0.05, 0.95, selection_str, 'Units', 'normalized', ...
        'FontSize', 16, 'Interpreter', 'none', 'Color', 'red', 'HorizontalAlignment', 'left');
    
    % Since sgtitle does not return a handle, set prompt_text to empty
    prompt_text = [];
end

function keyPressHandler(hObject, event)
    % Nested key press handler function

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
    selection_text = getappdata(hObject, 'selection_text');  % Retrieve selection text handle

    % Handle key presses
    if strcmp(event.Key, 'uparrow')
        disp('Up arrow pressed');
        arrow_key_input = 'up';
        gain = gain - 10;
        setappdata(hObject, 'gain', gain);
        hold off;
        % Re-plot with updated gain
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        % Retrieve current user_input for displaying
        current_input = getappdata(hObject, 'user_input');
        % Re-plot EEG with updated gain and existing user_input
        [~, selection_text_new] = plot_eeg_with_text(plot_values, fs, plot_labels, gain, current_input);

        % Update appdata with new text handles
        setappdata(hObject, 'prompt_text', []); % sgtitle doesn't return a handle
        setappdata(hObject, 'selection_text', selection_text_new);

    elseif strcmp(event.Key, 'downarrow')
        disp('Down arrow pressed');
        arrow_key_input = 'down';
        gain = gain + 10;
        setappdata(hObject, 'gain', gain);
        hold off;
        % Re-plot with updated gain
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        % Retrieve current user_input for displaying
        current_input = getappdata(hObject, 'user_input');
        % Re-plot EEG with updated gain and existing user_input
        [~, selection_text_new] = plot_eeg_with_text(plot_values, fs, plot_labels, gain, current_input);

        % Update appdata with new text handles
        setappdata(hObject, 'prompt_text', []); % sgtitle doesn't return a handle
        setappdata(hObject, 'selection_text', selection_text_new);

    elseif strcmp(event.Key, 'b') || strcmp(event.Key, 'c')
        disp(['Montage change requested: ', event.Key]);
        if strcmp(event.Key, 'b')
            montage_type = 'bipolar';
        elseif strcmp(event.Key, 'c')
            montage_type = 'car';
        end
        setappdata(hObject, 'montage_type', montage_type);
        hold off;
        % Re-plot with updated montage
        switch montage_type
            case 'bipolar'
                plot_values = bipolar_values;
                plot_labels = bipolar_labels;
            case 'car'
                plot_values = car_values;
                plot_labels = car_labels;
        end
        % Retrieve current user_input for displaying
        current_input = getappdata(hObject, 'user_input');
        % Re-plot EEG with updated montage and existing user_input
        [~, selection_text_new] = plot_eeg_with_text(plot_values, fs, plot_labels, gain, current_input);

        % Update appdata with new text handles
        setappdata(hObject, 'prompt_text', []); % sgtitle doesn't return a handle
        setappdata(hObject, 'selection_text', selection_text_new);

    elseif strcmp(event.Key, 'return')
        % User pressed Enter to confirm input
        if ismember(user_input, {'l', 'n', 'r', 'm'})
            disp(['User input confirmed: ', user_input]);
            % Store the user input
            setappdata(hObject, 'user_input', user_input);
            % Resume execution (but do NOT close the figure here)
            uiresume(hObject);
        else
            disp('Please press "l" (left), "r" (right), "m" (midline), or "n" (no spike) before pressing Enter to confirm.');
        end
    elseif ismember(event.Key, {'l', 'n', 'r', 'm'})
        % Store the user's input character but wait for Enter to confirm
        user_input = event.Key;
        disp(['You pressed "', user_input, '". Press Enter to confirm.']);
        % Store the tentative user input
        setappdata(hObject, 'user_input', user_input);
        
        % Update the selection display on the plot (Modified)
        set(selection_text, 'String', ['Selected: ', user_input]);
    else
        disp(['Other key pressed: ', event.Key]);
    end

    % Store the updated data back into appdata
    setappdata(hObject, 'arrow_key_input', arrow_key_input);
end

function [bipolar_values, bipolar_labels] = scalp_bipolar_clean(values, chLabels)
    % This function creates a bipolar montage
    
    montage = {'Fp1','F7';...
        'F7','T3';...
        'T3','T5';...
        'T5','O1';...
        '','';...
        'Fp2','F8';...
        'F8','T4';...
        'T4','T6';...
        'T6','O2';...
        '','';...
        'Fp1','F3';...
        'F3','C3';...
        'C3','P3';...
        'P3','O1';...
        '','';...
        'Fp2','F4';...
        'F4','C4';...
        'C4','P4';...
        'P4','O2';...
        '','';...
        'Fz','Cz';...
        '','';...
        '','';...
        'EKG1',''};
    
    bipolar_labels = cellfun(@(x,y) sprintf('%s-%s',x,y), montage(:,1), montage(:,2), 'UniformOutput', false);
    
    nsamples = size(values,1);
    nbi_channels = size(montage,1);
    bipolar_values = nan(nsamples, nbi_channels);
    
    % Convert ECGL and ECGR to EKG1 and EKG2
    if sum(strcmp(chLabels, 'ECGL')) ~= 0
        chLabels(strcmp(chLabels, 'ECGL')) = {'EKG1'};
    end
    if sum(strcmp(chLabels, 'ECGR')) ~= 0
        chLabels(strcmp(chLabels, 'ECGR')) = {'EKG2'};
    end
    
    for ib = 1:nbi_channels
        
        curr_mon = montage(ib,:);
        
        % Leave the empty rows as NaNs (just to designate blank space)
        if strcmp(curr_mon{1}, '') && strcmp(curr_mon{2}, '')
            continue;
        end
        
        % Get the channel labels that match
        lab1 = strcmp(chLabels, curr_mon{1});
        lab2 = strcmp(chLabels, curr_mon{2});
    
        if ~(sum(lab1) == 1 && sum(lab2) == 1)
            if contains(curr_mon{1}, 'EKG')
                if sum(lab1) == 0
                    continue;
                end
            else
                error('Missing expected channel: %s or %s', curr_mon{1}, curr_mon{2});
            end
        end
    
        if sum(lab2) == 0 && strcmp(chLabels(lab1), 'EKG1')
            bipolar_values(:, ib) = values(:, lab1);
            continue;
        end
        
        % Get difference in EEG signal
        val1 = values(:, lab1);
        val2 = values(:, lab2);
        bipolar_values(:, ib) = val1 - val2;
        
    end
end

function [car_values, car_labels] = car_montage_clean(values, chLabels)
    % This function creates a CAR (Common Average Reference) montage
    
    car_labels = chLabels;
    brain_channels = {'Fp1',...
        'F7',...
        'T3',...
        'T5',...
        'O1',...
        'Fp2',...
        'F8',...
        'T4',...
        'T6',...
        'O2',...
        'F3',...
        'C3',...
        'P3',...
        'F4',...
        'C4',...
        'P4',...
        'Fz',...
        'Cz'};
    
    % Take average of the brain channels
    is_brain = ismember(chLabels, brain_channels);
    avg = nanmean(values(:, is_brain), 2);
    
    % Do CAR (just average the non-skip channels)
    car_values = values; % Initialize
    car_values(:, is_brain) = values(:, is_brain) - repmat(avg, 1, sum(is_brain));
    
    for i = 1:length(chLabels)
        if is_brain(i) == 1
            car_labels{i} = [chLabels{i}, '-CAR'];
        end
    end
end

function plot_scalp_eeg_clean(values, fs, labels, gain)
    % This function plots the EEG data in a scalp montage
    
    added_offset = gain;
    
    dur = size(values,1)/fs;
    nchs = length(labels);
    
    offset = 0;
    
    ch_offsets = zeros(nchs,1);
    ch_bl = zeros(nchs,1);
    
    for ich = 1:nchs
        plot(linspace(0, dur, size(values,1)), values(:, ich) - offset, 'k');
        hold on;
        ch_offsets(ich) = offset;
        ch_bl(ich) = -offset + nanmedian(values(:, ich));
        offset = offset + added_offset;
        text(dur + 0.05, ch_bl(ich), sprintf('%s', labels{ich}), 'FontSize', 15);
    end
    plot([dur/2 dur/2], get(gca, 'ylim'), 'r--');
    
    yticklabels([]);
    xticks(1:floor(dur));
    ylim([-offset - added_offset, added_offset]);
    xlabel('Time (seconds)');
    set(gca, 'FontSize', 15);
    
    %% Add second markers
    for i = 1:floor(dur)
        plot([i, i], get(gca, 'ylim'), 'k--');
    end
end
