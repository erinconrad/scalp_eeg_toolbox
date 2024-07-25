function hdr = make_header_edf(values,labels,fs)

nsamples = size(values,1);
nchs = size(values,2);

hdr = edfheader("EDF+");
hdr.NumDataRecords = 1;
hdr.DataRecordDuration = seconds(nsamples/fs);
hdr.NumSignals = nchs;
hdr.SignalLabels = labels;
hdr.PhysicalDimensions = repelem("uV",nchs);
hdr.PhysicalMin = floor(min(values));
hdr.PhysicalMax = ceil(max(values));

hdr.PhysicalMin(isnan(hdr.PhysicalMin)) = -1e3;
hdr.PhysicalMax(isnan(hdr.PhysicalMax)) = 1e3;
    

hdr.DigitalMin = repmat(-32768,1,nchs);
hdr.DigitalMax = repmat(32767,1,nchs);



end