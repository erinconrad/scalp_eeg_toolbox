function values = bandpass_filter(values,fs)


d = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',70, ...
    'SampleRate',fs);

for i = 1:size(values,2)
    eeg = values(:,i);
    
    if sum(~isnan(eeg)) == 0
        continue
    end
    
    eeg(isnan(eeg)) = nanmean(eeg);
    eeg = filtfilt(d,eeg);   
    values(:,i) = eeg;
end


end