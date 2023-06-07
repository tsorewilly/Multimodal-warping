function [em, emg, force, glove] = staticAlignMultimodalData(dataIn)%em, emg, force, glove)
    %dataIn = [em, emg, force, glove];
    for i = 1 : size(dataIn, 2)
        CreateHands = double(split(string(dataIn(i).Creation), '.'));
        CreateTime(i) = (CreateHands(1)*3600 + CreateHands(2)*60 + CreateHands(3) + (CreateHands(4)/1000));        
        
        WriteHands = double(split(string(dataIn(i).Write), '.'));
        WriteTime(i) = (WriteHands(1)*3600 + WriteHands(2)*60 + WriteHands(3) + (WriteHands(4)/1000));
    end    
    
    [EarlyTime, EarlyIndex] = sort(CreateTime, 'descend');
    LeadTime = abs(EarlyTime - EarlyTime(size(EarlyTime, 2)));
    
    [LateTime, LateIndex] = sort(WriteTime, 'ascend');
    LagTime = abs(LateTime - LateTime(1));
    
    for i = 1 : size(EarlyTime, 2)
        recRows = size(dataIn(EarlyIndex(i)).Data,1);
        recPeriod = dataIn(EarlyIndex(i)).Diff;
        
        dataIn(EarlyIndex(i)).StartTime = CreateTime(EarlyIndex(i));
        dataIn(EarlyIndex(i)).EndTime = WriteTime(EarlyIndex(i));
        dataIn(EarlyIndex(i)).LeadTime = LeadTime(EarlyIndex(i));
        dataIn(EarlyIndex(i)).LagTime = LagTime(LateIndex(EarlyIndex(i)));
                
        %dataIn(EarlyIndex(i)).AlignStartOld = floor((recRows * LeadTime(EarlyIndex(i)))/recPeriod);        
        dataIn(EarlyIndex(i)).AlignStart = floor((recRows * (CreateTime(EarlyIndex(1)) - CreateTime(EarlyIndex(i))))/recPeriod)+1; %Added 1 to avoid zero indexing        
        dataIn(EarlyIndex(i)).AlignEnd = floor((recRows * (dataIn(EarlyIndex(i)).Diff - (WriteTime(EarlyIndex(i)) - WriteTime(LateIndex(1)))))/recPeriod);
        %dataIn(EarlyIndex(i)).AlignEndOld = floor((recRows * dataIn(EarlyIndex(1)).Diff - LeadTime(1))/recPeriod);
    end 
   
    em = dataIn(1);
    emg = dataIn(2);
    force = dataIn(3);
    glove = dataIn(4);
end

