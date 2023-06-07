function VCF = vcf1(signal,fs)
    t = abs(fft(signal));
    t(1) = [];
    tmax = max(t);
    t(ceil(end/2):end) = [];
    abovecutoff = t > tmax / 2;   %3 dB is factor of 2
    lowbin  = find(abovecutoff, 1, 'first');
    highbin = sum(abovecutoff);
    centbin = sqrt(lowbin * highbin);   %geometric mean
    cenFreqBin = floor(centbin);
end
%Now you know that floor(centbin) is the fft bin containing the center frequency, but you don't know the center frequency itself. 
%To know the frequency itself you need to convert between bin numbers and frequencies,
%which depends upon your sample resolution -- which is something that cannot be know from just the signal itself.

%Here I am using "fft bin number" in the sense of bin #1 corresponding to 1 cycle/period, bin #2 corresponding to 2 cycles / period, and so on.
%Having code on hand you should now be able to refer back to your textbook 
%check definitions of fft and of centre frequencies and cutoff frequencies in order to figure out how the code works. 