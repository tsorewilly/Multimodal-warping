function dataOut = updateData(em, emg, force, glove)
    %emLength = size(em, 1);
    emNorm = em.Diff/size(em.Data, 1);   %Normalize the em data using its duration vs length   
    %%
    %Convert emg.Creation time into seconds and add glove signal duration
    hands = double(split(string(glove.Creation), '.'));
    glove.Write = datetime((hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000)) + glove.Diff, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    %%
    %convert emgLength into appropriate duration using its resolution frequency
    emgAt40Hz = downsample(emg.Data, floor(emg.Resolution/em.Resolution));
    emg.Diff = size(emgAt40Hz, 1) * emNorm;
    
    %Convert emg.Write time into seconds and add the signal duration
    hands = (double(split(string(em.Creation), '.')) + double(split(string(glove.Creation), '.')))/2;
    emg.Creation = datetime((hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000)), 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    emg.Write = datetime((hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000)) + emg.Diff, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    %%
    %convert emgLength into appropriate duration using its resolution frequency
    force.Diff = ((size(force.Data, 1) * em.Resolution)/force.Resolution) * emNorm;
    
    %Convert emg.Write time into seconds and add the signal duration
    hands = (double(split(string(emg.Creation), '.')) + double(split(string(glove.Creation), '.')))/2;
    force.Creation = datetime((hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000)), 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    force.Write = datetime((hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000)) + force.Diff, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    
    %%
    %Add a similar timezone to files that doesn't have it.
    em.Creation = datetime(em.Creation, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    em.Write = datetime(em.Write, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    glove.Creation = datetime(glove.Creation, 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
    dataOut = [em, emg, force, glove];
end