% validate spike detections
%{
Instructions:
1. Install dependencies (anything??) and add to path
2. Update paths below
3. Navigate to directory containing this script and type
  >> clean_spike_validator

The script will loop over all edf files in a single directory and display
the EEG in each edf file. The user will be prompted to select if the spike
is left (l), right (r), midline (m), or if there is no spike (n). The user
will then press return/enter to confirm the choice. Selections will be
saved in a csv file in the directory containing edf files.

The montage can be toggled between bipolar (b) and CAR (c).

Pressing the up or down array will adjust the gain.

The default overwrite setting (0) is to not re-present any edf files the
user has already completed. Changing this to overwrite = 1 will allow the
user to repeat spike selections for all files.

%}

clear

%% Paths - edit this!!!
edf_dir = '/Users/erinconrad/Desktop/research/scalp_eeg_toolbox/results/example_edfs/';
out_file = 'selections.csv';

%% overwrite settings - edit this as needed
% overwrite = 0 (default): skip previously completed files
% overwrite = 1: re-write selections of all files
% overwrite = 2: show all files, but do not save any selections
overwrite = 0;


%{ ----- NO NEED TO EDIT BELOW THIS ------------ %}


%% Define variable names in table
% Define the variable names
variableNames = {'Filename', 'Response'};  % Replace with your desired variable names

%% Load output file if it exists and initialize table
if overwrite == 0
    % Load the out file if it exists
    if exist(fullfile(edf_dir,out_file),"file") ~= 0
        responses = readtable(fullfile(edf_dir,out_file),...
            'ReadVariableNames', true,'HeaderLines', 0);
    else
        responses = table('Size', [0, 2], 'VariableTypes', {'string', 'string'}, ...
                  'VariableNames', variableNames);
    end
else
    responses = table('Size', [0, 2], 'VariableTypes', {'string', 'string'}, ...
                  'VariableNames', variableNames);
end


% Get the list of edf files
edf_files = dir(fullfile(edf_dir, '*.edf'));

% Loop through each EDF file
for i = 1:length(edf_files)

    %% Figure out the file path and whether to run it
    % Get the full path of the current EDF file
    edf_filename = fullfile(edf_dir, edf_files(i).name);

    % Look and see if we've done it already (if not overwriting)
    if overwrite == 0 || overwrite == 2
        completed_files = responses.Filename;

        if any(strcmp(completed_files,edf_files(i).name))
            fprintf('\nAlready did %s, skipping\n',...
                edf_files(i).name);
            continue;
        end
    end

    
    
    %% Load and process the EDF file
    [values, chLabels, fs] = read_in_edf(edf_filename);

    % Initialize gain
    gain = 30;
    
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
    [bipolar_values, bipolar_labels] = scalp_bipolar_clean(values,chLabels);
   
    % also get car
    [car_values,car_labels] = car_montage_clean(values,chLabels);
    
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
    figure
    hFig = gcf;
    set(hFig, 'Position', [10 10 1400 1000]);
    tiledlayout(1, 1, 'Padding', 'compact');
    nexttile
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
    fprintf('\nDisplaying %s\n',(edf_files(i).name))
    
    %% Handle user input
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
    
    %% Handle end of file and write responses
    % Now you can safely close the figure
    if ishandle(hFig)
        close(hFig);
    end
    
    % Store the filename and response in the responses array
    responses(end+1,:) = {edf_files(i).name, user_input};
    
    if overwrite ~=2
        % Save the table as a CSV file (save as you go)
        writetable(responses, fullfile(edf_dir, out_file));
    end

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

bipolar_labels = cellfun(@(x,y) sprintf('%s-%s',x,y),montage(:,1),montage(:,2),'uniformoutput',false);

nsamples = size(values,1);
nbi_channels = size(montage,1);
bipolar_values = nan(nsamples,nbi_channels);

% convert EKGL and R to 1 and 2
if sum(strcmp(chLabels,'ECGL')) ~= 0, chLabels(strcmp(chLabels,'ECGL')) = {'EKG1'}; end
if sum(strcmp(chLabels,'ECGR')) ~= 0, chLabels(strcmp(chLabels,'ECGR')) = {'EKG2'}; end

for ib = 1:nbi_channels
    
    curr_mon = montage(ib,:);
    
    % leave the empty rows nans (just to designate blank space)
    if strcmp(curr_mon{1},'') && strcmp(curr_mon{2},'')
        continue;
    end
    
    % get the channel labels that match
    lab1 = strcmp(chLabels,curr_mon{1});
    lab2 = strcmp(chLabels,curr_mon{2});

    
    if ~(sum(lab1)==1 && sum(lab2)==1)
        if contains(curr_mon{1},'EKG')
            if sum(lab1) == 0
                continue
            end
        else
            error('missing expected channel')
        end
    end

    if sum(lab2) == 0 && strcmp(chLabels(lab1),'EKG1')
        bipolar_values(:,ib) = values(:,lab1);
        continue
    end
    
    % get difference in eeg signal
    val1 = values(:,lab1);
    val2 = values(:,lab2);
    bipolar_values(:,ib) = val1-val2;
    
    
end


end


function [car_values,car_labels] = car_montage_clean(values,chLabels)
% This function creates a car montage

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
is_brain = ismember(chLabels,brain_channels);
avg = nanmean(values(:,is_brain),2);

% Do car (just average the non-skip chs)
car_values(:,is_brain) = values(:,is_brain) - repmat(avg,1,sum(is_brain));

for i = 1:length(chLabels)
    if is_brain(i) == 1
        car_labels{i} = [chLabels{i},'-CAR'];
    end
end

end

function plot_scalp_eeg_clean(values,fs,labels,gain)
% This function plots the eeg data

added_offset = gain;

dur = size(values,1)/fs;
nchs = length(labels);

offset = 0;

ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for ich = 1:nchs
    plot(linspace(0,dur,size(values,1)),values(:,ich) - offset,'k');
    hold on
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    offset = offset+added_offset;
    text(dur+0.05,ch_bl(ich),sprintf('%s',labels{ich}),'fontsize',15)
    
    %{
    if ich < nchs
        if ~isnan(min(values(:,ich)) - max(values(:,ich+1)))
            offset = offset - (min(values(:,ich)) - max(values(:,ich+1)));
        end
        
    end
    %}
    
end
plot([dur/2 dur/2],get(gca,'ylim'),'r--')

yticklabels([])
xticks(1:floor(dur))
ylim([-offset-added_offset added_offset])
xlabel('Time (seconds)')
set(gca,'fontsize',15)

%% Add second markers
for i = 1:dur
    plot([i,i],get(gca,'ylim'),'k--');
    
end


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
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
    
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
    plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
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
        plot_scalp_eeg_clean(plot_values, fs, plot_labels, gain);
elseif strcmp(event.Key, 'return')
    % User pressed Enter to confirm input
    if strcmp(user_input, 'l') || strcmp(user_input, 'n') || strcmp(user_input, 'r') || strcmp(user_input, 'm')
        disp(['User input confirmed: ', user_input]);
        % Store the user input
        setappdata(hObject, 'user_input', user_input);
        % Resume execution (but do NOT close the figure here)
        uiresume(hObject);
    else
        disp('Please press "l" (left), "r" (right), "m" (midline), or "n" (no spike) before pressing Enter to confirm.');
    end
elseif strcmp(event.Key, 'l') || strcmp(event.Key, 'n') || strcmp(event.Key, 'r') || strcmp(event.Key, 'm')
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
