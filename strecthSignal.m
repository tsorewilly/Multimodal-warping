function StretchedSignal = strecthSignal(signal_cells, signal_id, ref_id, shift_method, w)
    if(~exist('shift_method','var'))
        shift_method = 0; %No Need to use a particular shift method
    end
    if(~exist('w','var'))
        w = 0.5; %Use default warping window
    end
    %Divide signal into two parts wrt ref_signal
    ref_len = length(signal_cells{ref_id});
    signal = signal_cells{signal_id};
    shift_len = length(signal);
    
    warpLength = ref_len - shift_len;
    fillIns = (randperm(ref_len, ref_len - shift_len))';                    %Select random indices to stretch signal    
    
    if shift_method == 1
        num = ((abs(fillIns - ceil(ref_len/2)))./w);
        fillIns = num./((1+num).^2);
    elseif  shift_method == 2
        fillIns = 1./(abs(fillIns - ceil(ref_len/2)).^2);                   %Used in Keogh 2001   
    end

    fillIns = sort(fillIns);
    
    for j = 1 : warpLength
        signal(fillIns(j)+1:end+1, :) = signal(fillIns(j):end, :);          %Shift the leading signal forward        
        %signal(fillIns(j), :) = zeros(1,size(signal,2));                   %Replace value at the index with Zero
        if j == 1
            dsignal = (signal(1 : fillIns(j), :));
        elseif j == warpLength
            dsignal = (signal(fillIns(j-1) : fillIns(j), :));
        else
            dsignal = (signal(fillIns(j) : fillIns(j+1), :));
        end
        signal(fillIns(j), :) = mean(dsignal) + std(dsignal);
        %signal_Left = signal(1:floor(Ref_Length/2));
        %signal_Right = signal(floor(Ref_Length/2)+1:end);
        %Strecth Signals with Right Alignment 
        %Strecth Signals with Left Alignment 
    
    
    end
    StretchedSignal = signal;
end
