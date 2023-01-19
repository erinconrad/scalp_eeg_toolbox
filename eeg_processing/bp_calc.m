function all_power = bp_calc(run_values,fs,skip)

freqs = [1 4;...
    4 8;...
    8 12;...
    12 30;...
    30 70];
nfreqs = size(freqs,1);

nchs = size(run_values,2);
all_power = nan(nchs,nfreqs);

% turn skip channels into zeros
more_skip = ismember(1:nchs,skip) | all(isnan(run_values),1);
run_values(:,more_skip) = 0;

for f = 1:nfreqs
    curr_freq = freqs(f,:);
    p = bandpower(run_values,fs,curr_freq);
    all_power(:,f) = p;

end

all_power(more_skip,:) = nan;


end