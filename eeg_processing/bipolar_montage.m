function [values,clean_labels,bipolar_labels,chs_in_bipolar,which_chs_bipolar,mid_locs,mid_anatomy] =...
    bipolar_montage(values,chLabels,which_chs,locs,anatomy,name)

%{
This function takes a chunk of multi-channel EEG data, along with channel
names, and outputs the data in a bipolar montage. The output values(:,ich)
equals the input values(:,ich) - values(:,jch), where jch is the next
numbered contact on the same electrode as ich. For example. if ich is RA1,
then jch is RA2 and values(:,RA1) will be the old values(:,RA1) - old
values(:,RA2). If ich is the last contact on the electrode, then
values(:,ich) is defined to be nans.
%}

%% Initialize output variables
nchs = size(values,2);
chs_in_bipolar = nan(nchs,2);
old_values = values;
which_chs_bipolar = [];

%% Input variable handling
if isempty(which_chs)
    which_chs = 1:size(values,2);
end

%% Decompose chLabels
[clean_labels,elecs,numbers] = decompose_labels(chLabels,name);
bipolar_labels = cell(nchs,1);

%% Bipolar montage
for ch = 1:nchs
    
    % Initialize it as nans
    out = nan(size(values,1),1);
    
    % Get the clean label
    label = clean_labels{ch};

    % get the non numerical portion of the electrode contact
    label_non_num = elecs{ch};

    % get numerical portion
    label_num = numbers(ch);

    % see if there exists one higher
    label_num_higher = label_num + 1;
    higher_label = [label_non_num,sprintf('%d',label_num_higher)];
    if sum(strcmp(clean_labels(:,1),higher_label)) > 0
        higher_ch = find(strcmp(clean_labels(:,1),higher_label));
        out = old_values(:,ch)-old_values(:,higher_ch);
        bipolar_label = [label,'-',higher_label];
        chs_in_bipolar(ch,:) = [ch,higher_ch];
        
        % Only make this a run channel if both members of the bipolar
        % channel were run channels
        if ismember(ch,which_chs) && ismember(higher_ch,which_chs)
            which_chs_bipolar = [which_chs_bipolar;ch];
        end
        
    elseif strcmp(label,'FZ') % exception for FZ and CZ
        if sum(strcmp(clean_labels(:,1),'CZ')) > 0
            higher_ch = find(strcmp(clean_labels(:,1),'CZ'));
            out = old_values(:,ch)-old_values(:,higher_ch);
            bipolar_label = [label,'-','CZ'];
            chs_in_bipolar(ch,:) = [ch,higher_ch];
            
            if ismember(ch,which_chs) && ismember(higher_ch,which_chs)
                which_chs_bipolar = [which_chs_bipolar;ch];
            end
        end
        
    else
        % allow it to remain nans
        bipolar_label = '-';
    end
    values(:,ch) = out;
    bipolar_labels{ch} = bipolar_label;
    
    
    
end

which_chs_bipolar = unique(which_chs_bipolar);

%% Get location of midpoint between the bipolar channels
if ~isempty(locs)
    mid_locs = nan(length(bipolar_labels),3);
    mid_anatomy = cell(length(bipolar_labels),1);
    for i = 1:length(bipolar_labels)

        % Get the pair
        ch1 = chs_in_bipolar(i,1);
        ch2 = chs_in_bipolar(i,2);

        if isnan(ch1) || isnan(ch2)
            continue
        end

        % get the locs
        loc1 = locs(ch1,:);
        loc2 = locs(ch2,:);

        % get midpoint
        midpoint = (loc1 + loc2)/2;
        mid_locs(i,:) = midpoint;


    end
else
    mid_locs = [];
end

%% Get anatomy
if ~isempty(anatomy)
    mid_anatomy = cell(length(bipolar_labels),1);
    for i = 1:length(bipolar_labels)

        % Get the pair
        ch1 = chs_in_bipolar(i,1);
        ch2 = chs_in_bipolar(i,2);

        if isnan(ch1) || isnan(ch2)
            continue
        end

        
        % get anatomy of each
        anat1 = anatomy{ch1};
        anat2 = anatomy{ch2};
        midanat = [anat1,'-',anat2];
        mid_anatomy{i} = midanat;

    end
else
    mid_anatomy = [];
end


end