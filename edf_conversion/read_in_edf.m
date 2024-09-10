function [out,chLabels,fs] = read_in_edf(filepath)

data = edfread(filepath);
info = edfinfo(filepath);
fs = info.NumSamples(1);

chs = info.SignalLabels;
nchs = length(chs);

nsamples = size(data,1) * fs;

out = nan(nsamples,nchs);

for i = 1:nchs
    out(:,i) = cell2mat(data.(chs(i)));

end

chLabels = cellstr(chs);

end