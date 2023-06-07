function mdf1(Signal, Fs)
    Fs = 1024;
    x = Signal;%randn(841700,1);
    window = 100;
    [S,F,T,P] = spectrogram(x,window,window/2,window,Fs);
    for nn = 1:size(P,2)
    normcumsumpsd = cumsum(P(:,nn))./sum(P(:,nn));
    Ind = find(normcumsumpsd <=0.5,1,'last');
    medianfreqs(nn) = F(Ind);
    end
    plot(T,medianfreqs);
    xlabel('Time (seconds)'); 
    ylabel('Median Frequency (Hz)');