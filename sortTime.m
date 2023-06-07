function data = sortTime(data, opMode, sortMode)
    %data: array of start or stop times for a set og signals to be aligned
    %min = 1; max = 2; mean = 3; sort = 4;
    if(~exist('sortMode','var')); sortMode='descend';   end
    for i = 1 : size(data, 2)
        hands = double(split(string(data(i).Creation), '.'));
        dsTime(i) = (hands(1)*3600 + hands(2)*60 + hands(3) + (hands(4)/1000));
    end
    
    if (opMode == 1)
        EarlyTime = min(dsTime);
    elseif (opMode == 2)
        EarlyTime = max(dsTime);
    elseif (opMode == 3)
        EarlyTime = mean(dsTime);
    elseif (opMode == 4)
        [EarlyTime, TimeIndex] = sort(dsTime, sortMode);
        lagTime = abs(EarlyTime - EarlyTime(1));
    end
    
    for i = 1 : size(EarlyTime, 2)
        %sortedTime(i) = datetime(EarlyTime(i), 'ConvertFrom', 'posixtime', 'TimeZone','Etc/GMT','TicksPerSecond',1000, 'Format','HH.mm.ss.SSS');
        data(TimeIndex(i)).LagTime = lagTime(i);
    end 
end