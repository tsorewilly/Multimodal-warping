function fileTimes = getFileTimes(filename)
    fileTimes = GetFileTime(filename, 'Windows');
    %startTime = fileTimes.Write;
    fileTimes.Creation = datetime(fileTimes.Creation,'Format','HH.mm.ss.SSS');      
    hands = double(split(string(fileTimes.Creation), '.'));
    tm1 = hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000);

    %endTime = fileTimes.Write;
    fileTimes.Write = datetime(fileTimes.Write,'Format','HH.mm.ss.SSS');  
    hands = double(split(string(fileTimes.Write), '.'));
    tm2 = hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000);
    
    diff = tm2 - tm1;
    
    fileTimes.Diff = abs(diff);
    
    fileTimes = rmfield(fileTimes,'Access');
end