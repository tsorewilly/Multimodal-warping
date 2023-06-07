function s = DistShiftSeq(sequence, w, method)
    if method == 1
        num = exp(-(abs(sequence(:) - sequence(ceil(size(sequence,2)/2))))./w);
        s = num./((1+num).^2);
    elseif  method == 2
        s = 1./(abs(sequence(:) - sequence(ceil(size(sequence,2)/2))).^2); %Used in Keogh 2001
end