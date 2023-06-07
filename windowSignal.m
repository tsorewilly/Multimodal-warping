function [TRs, VRs] = windowSignal(signal, n, res)
        
     RefSigLen = size(signal,1); 
     TWC = 0; %TotalwindowCount
     TRs = zeros(1, floor(RefSigLen/n) + 1);
     for j = 1 : n        
         %window = n*j-n+1 : n*j;
         %((RefSigLen - TWC)/(n-j))
         %RefSigLen - TWC
         %n-j
         %if mod(((RefSigLen - TWC)/(n-j)), 2) > 0
         if ((RefSigLen - TWC)/(n-j)) > floor(RefSigLen/n)
            lent = floor(RefSigLen/n) + 1;
         else
            lent = floor(RefSigLen/n);            
         end
         TWC = TWC + lent;
         if TWC > RefSigLen
            %TRs = zeros(1, floor(RefSigLen/n) + 1);
            TRs(j,1:length([TRs(end)+1 : RefSigLen])) = [TRs(end)+1 : RefSigLen];
         else
            TRs(j,:) = [TRs(end)+1 : TWC];
         end
         
     end
     %Fill the columns with zero index on last row with the last index 
     TRs(TRs==0) = max(max(TRs));
     %movingObject(movingObject<=-40000)=1;  
%     Get the corresponding amplitudes
     for j = 1 : n
         VRs(j,:) = (signal(TRs(j,:)))';
         %VRs(j,:) = signal([TRs(j, 1) TRs(j, end)]);
     end
     
%     sig_per_act = n * res;  
% 
%     RefSigLen = size(signal,1);
%     
%     for j = 1 : RefSigLen/sig_per_act        
%         window = sig_per_act*j-sig_per_act+1 : sig_per_act*j;
%         TRs(j,:) = [window(1) window(end)];
%     end
% 
%     TRs(j+1,:) = [window(end)+1 length(signal)]; %Fill up the remaining signal
% 
%     Get the corresponding amplitudes
%     for j = 1 : length(TRs)
%         if diff(TRs(j, 1:2)) < sig_per_act-1
%             VRs(j,1 : mod(RefSigLen, sig_per_act)) = signal(TRs(j, 1) : TRs(j, 2));
%             VR_times(j,1 : mod(RefSigLen, sig_per_act)) = [TRs(j, 1) : TRs(j, 2)];
%         else
%             VRs(j,:) = signal(TRs(j, 1) : TRs(j, 2));
%             VR_times(j,:) = [TRs(j, 1) : TRs(j, 2)];
%         end
%     end
end
