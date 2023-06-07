function MPF = mpf1(signal,fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
	% Average and take the square root of each window
	%y(index) = signal(i:i+windowlength-1);
    %average each window
    %y(index) = mean(signal(i:i+windowlength-1));
    L1 = length(signal);
    cx1 = xcorr(signal,'unbiased');
    cxk1 = fft(cx1,L1);
    px1 = abs(cxk1);
    f1 = (0:L1-1)*fs/L1;
    df1 = fs/L1;
    p1=(sum(px1(1:L1/2-1))+sum(px1(1:L1/2)))/2.*df1;
    pf1=(sum(px1(1:L1/2-1).*[1:L1/2-1]'.*df1)...
            + sum(px1(1:L1/2).*[1:L1/2]'.*df1)...
        )/2*df1;
    MPF = pf1/p1;
end

