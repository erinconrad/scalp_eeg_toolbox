%{
This code takes in an edf file with the path "filepath" and outputs the
time series data, channel labels, and sampling rate.

Inputs:
- filepath: the path to the edf file

Outputs:
out: an nsamples x nchannels matrix containing time series data
chLabels: a nchannels x 1 cell containing channel labels
fs: the sampling frequency
%}

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