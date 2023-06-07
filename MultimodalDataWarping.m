%vars2keep = {'Cnt', 'SpectSegMusc'};
%clearvars('-except',vars2keep{:});
clc; clear;
action = "PUSH";
if(~exist('Cnt','var'))
    Cnt = 0; %Use default warping window
end
%Abductor pollicis brevis (APB), Flexor carpi radialis (FCR), Dorsal interossei (DI), and Extensor carpi radialis (ECR)
frcFolder = 'E:\Rabbit\2020.05.28\omisore-02\force data\';
emgFolder = 'E:\Rabbit\2020.05.28\omisore-02\emg data\';
emFolder =  'E:\Rabbit\2020.05.28\omisore-02\EM data\';
glvFolder = 'E:\Rabbit\2020.05.28\omisore-02\glove data\';

%% 1) Load Timeseries signals
%Load Raw sEMG signal
fileEMG = uigetfile('*.acq', 'Select Signal To Process', emgFolder);

if(fileEMG~=0)
    emg = getFileTimes([emgFolder, fileEMG]);
    emg.Resolution = 1024;
    emgFile = load_acq([emgFolder, fileEMG]);
    emg.Data = emgFile.data;
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('Select A File of EMG Signal','File Selection Error',opts);
    error('Select A File of EMG Signal');
end

%Load Raw Force Signal
fileForce = uigetfile('*.txt', 'Select Signal To Process', frcFolder);
if(fileForce~=0)
    force = getFileTimes([frcFolder, fileForce]);
    force.Resolution = 60;
    fFile = fopen([frcFolder, fileForce]);
    force.Data = textscan(fFile,'%s %s %f32 %s');
    force.Data =  force.Data{3};
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('Select A Recorded File for Force Data','File Selection Error',opts);
    error('Select A Recorded File for Force Data');
end

%Load EM Motion Signal
fileEM = uigetfile('*.csv', 'Select Signal To Process', emFolder);
if(fileEM~=0)
    em = getFileTimes([emFolder, fileEM]);
    em.Resolution = 40;
    opts = delimitedTextImportOptions("NumVariables", 25);    
    opts.DataLines = [2, Inf]; opts.Delimiter = ",";% Specify range and delimiter

    % Specify column names and types
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Q0", "Qx", "Qy", "Qz", "Tx", "Ty", "Tz", "Error", "Var14", "Var15", "Var16", "Var17", "Q1", "Qx1", "Qy1", "Qz1", "Tx1", "Ty1", "Tz1", "Error1"];
    opts.SelectedVariableNames = ["Q0", "Qx", "Qy", "Qz", "Tx", "Ty", "Tz", "Error", "Q1", "Qx1", "Qy1", "Qz1", "Tx1", "Ty1", "Tz1", "Error1"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 14, 15, 16, 17], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2, 3, 4, 5, 14, 15, 16, 17], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";    opts.EmptyLineRule = "read";

    em.Data = table2array(readtable([emFolder, fileEM], opts)); % Import the data    
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('Select A Recorded File for EM Data','File Selection Error',opts);
    error('Select A Recorded File for EM Data');
end

%Load Glove Motion Signal
fileGlv = uigetfile('*.csv', 'Select Signal To Process', glvFolder);
if(fileGlv~=0)
    glv = getFileTimes([glvFolder, fileGlv]);
    glv.Resolution = 55;
    opts = delimitedTextImportOptions("NumVariables", 30);
    opts.DataLines = [1, Inf];    opts.Delimiter = ","; % Specify range and delimiter    

    % Specify column names and types
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "ThumbNear", "ThumbFar", "ThumbIndex", "IndexNear", "IndexFar", "IndexMiddle", "MiddleNear", "MiddleFar", "MiddleRing", "RingNear", "RingFar", "RingLittle", "LittleNear", "LittleFar"];   
    opts.SelectedVariableNames = ["Var1", "ThumbNear", "ThumbFar", "ThumbIndex", "IndexNear", "IndexFar", "IndexMiddle", "MiddleNear", "MiddleFar", "MiddleRing", "RingNear", "RingFar", "RingLittle", "LittleNear", "LittleFar"];                    
    opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double","double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";    opts.EmptyLineRule = "read";

    glv.Data = table2array(readtable([glvFolder, fileGlv], opts));    % Import the whole glove data
    lent = size(glv.Data, 1); cond=1;   duration=0;
    while cond==1
       if isnan(double(glv.Data(lent,1)))
           opts.DataLines = [lent+1, Inf];
           gloveStart =  split(glv.Data(lent-3,1), 'at ');
           gloveStart = datetime(gloveStart(2),'Format','HH.mm.ss.SSS');
           cond=0;
       else
           duration = duration + double(glv.Data(lent,1));
           lent=lent-1;
       end
    end

    opts.SelectedVariableNames = ["ThumbNear", "ThumbFar", "ThumbIndex", "IndexNear", "IndexFar", "IndexMiddle", "MiddleNear", "MiddleFar", "MiddleRing", "RingNear", "RingFar", "RingLittle", "LittleNear", "LittleFar"];
    glv.Data = table2array(readtable([glvFolder, fileGlv], opts));    % Import just this session data
    glv.Creation = gloveStart;     glv.Diff = duration;
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('Select A Recorded File for Glove Data','File Selection Error',opts);
    error('Select A Recorded File for Glove Data');
end
clear opts;

%% 2) Warping: Reference signal Selection, Streching, and Filling

%Select Reference signal (R) and resample all data base on R's frequency 
Ref_resolution = emg.Resolution; %(Set reference as EM)
sampled_EMG = resample(emg.Data, Ref_resolution, emg.Resolution);
sampled_EM = resample(cellfun(@(x)str2double(x), em.Data), Ref_resolution, em.Resolution);
sampled_EM(sampled_EM < -200) = 0.01; sampled_EM(sampled_EM > 200) = 0.01; 
sampled_glv = resample(glv.Data, Ref_resolution, glv.Resolution);
sampled_force = resample(double(force.Data), Ref_resolution, force.Resolution);

figure(1)
subplot(4,4,1);plot(emg.Data);  subplot(4,4,2); plot(str2double(em.Data)); 
subplot(4,4,3);plot(glv.Data);  subplot(4,4,4); plot(force.Data);

subplot(4,4,5);plot(sampled_EMG(:,1));  subplot(4,4,6);plot(sampled_EM(:,4)); 
subplot(4,4,7);plot(sampled_glv(:,1));  subplot(4,4,8);plot(sampled_force(:,1));

%%
ALL_Signals = {sampled_EMG, sampled_EM, sampled_glv, sampled_force};
ALL_Signals_resolution = {emg.Resolution, em.Resolution, glv.Resolution, force.Resolution};
StretchedSignal = ALL_Signals;
[Ref_Length, Ref_Index] = max(cellfun(@length, ALL_Signals));

for k = 1:length(ALL_Signals)    
    if k ~= Ref_Index
        %StretchedSignal{k} = strecthSignal(ALL_Signals, k, Ref_signal);                        %Strecth Signals to Alignment
        shift_len = length(ALL_Signals{k});
        warpLength = Ref_Length - shift_len;
        fillIns = sort(randperm(Ref_Length, Ref_Length - shift_len))';                          %Select random indices to stretch signal    
        for dataId = 1 : warpLength
            StretchedSignal{k}(fillIns(dataId)+1:end+1, :) = StretchedSignal{k}(fillIns(dataId):end, :);  %Shift the leading signal forward                    
            if dataId <= 6 || warpLength-dataId <= 6
                rearAs = 6;
            else
                rearAs = dataId;
            end
            dsignal = (StretchedSignal{k}(fillIns(rearAs-5) : fillIns(rearAs+5), :));
            StretchedSignal{k}(fillIns(dataId), :) = mean(dsignal) + std(dsignal);
            %StretchedSignal{k}(fillIns(j), :) = zeros(1,size(signal,2));                       %Replace value at the index with Zero
        end
    end
end

subplot(4,4,9);plot(StretchedSignal{1}(:,1)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Abductor pollicis brevis EMG');
subplot(4,4,10);plot(StretchedSignal{1}(:,2)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Flexor carpi radialis EMG');
subplot(4,4,11);plot(StretchedSignal{1}(:,3)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Dorsal interossei EMG');
subplot(4,4,12);plot(StretchedSignal{1}(:,4)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Extensor carpi radialis EMG'); 

    [halfFFT1, offset_signal1, medfreq1] = process_emg(StretchedSignal{1}(:,1), emg.Resolution);  %ALL_Signals_resolution{Ref_Index}
        subplot(4,4,13);plot(halfFFT1(:,1), halfFFT1(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT2, offset_signal2, medfreq2] = process_emg(StretchedSignal{1}(:,2), emg.Resolution); 
        subplot(4,4,14);plot(halfFFT2(:,1), halfFFT2(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT3, offset_signal3, medfreq3] = process_emg(StretchedSignal{1}(:,3), emg.Resolution); 
        subplot(4,4,15);plot(halfFFT3(:,1), halfFFT3(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT4, offset_signal4, medfreq4] = process_emg(StretchedSignal{1}(:,4), emg.Resolution); 
        subplot(4,4,16);plot(halfFFT4(:,1), halfFFT4(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    %The apart above, signal offset removal and nfft, was not used.

    
%% 3) Time-series Processing for Feature Extraction 
%Select the pre-filtered sEMG

figure(2);
[FilteredEMGData(:,1), EMG_Segments(:,:,1), RMS(1,:), Xscale(1,:), rmsXscale(1,:)] = filtEMGSignal(StretchedSignal{1}(:,5), emg.Resolution, 50, 500, 64);    
    subplot(411); plot(Xscale(1,:),abs(StretchedSignal{1}(:,1)), '-b'); hold on; plot (rmsXscale(1,:), RMS(1,:), '-r', 'linewidth', 2); axis tight
    xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
[FilteredEMGData(:,2), EMG_Segments(:,:,2), RMS(2,:), Xscale(2,:), rmsXscale(2,:)] = filtEMGSignal(StretchedSignal{1}(:,6), emg.Resolution, 50, 500, 64);    
    subplot(412); plot(Xscale(2,:),abs(StretchedSignal{1}(:,2)), '-b'); hold on; plot (rmsXscale(2,:), RMS(2,:), '-r', 'linewidth', 2); axis tight
    xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
[FilteredEMGData(:,3), EMG_Segments(:,:,3), RMS(3,:), Xscale(3,:), rmsXscale(3,:)] = filtEMGSignal(StretchedSignal{1}(:,7), emg.Resolution, 50, 500, 64);    
    subplot(413); plot(Xscale(3,:),abs(StretchedSignal{1}(:,3)), '-b'); hold on; plot (rmsXscale(3,:), RMS(3,:), '-r', 'linewidth', 2); axis tight
    xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
[FilteredEMGData(:,4), EMG_Segments(:,:,4), RMS(4,:), Xscale(4,:), rmsXscale(4,:)] = filtEMGSignal(StretchedSignal{1}(:,8), emg.Resolution, 50, 500, 64);    
    subplot(414); plot(Xscale(4,:),abs(StretchedSignal{1}(:,4)), '-b'); hold on; plot (rmsXscale(4,:), RMS(4,:), '-r', 'linewidth', 2); axis tight
    xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);

%ThumbFinger = deNoiseEM(double(StretchedSignal{1,2}(:, 1:7)), 200);
%IndexFinger = deNoiseEM(double(StretchedSignal{1,2}(:, 9:15)), 200);
%interHand3D = distMoved(ThumbFinger(:,5:7), IndexFinger(:,5:7))';
%intraThumb3D = distMoved(ThumbFinger(:,5:7))';
%intraIndex3D = distMoved(IndexFinger(:,5:7))';

EM1_rot_mat = quat2rotm(StretchedSignal{1,2}(:,1:4));
EM1_trans_vec = StretchedSignal{1,2}(:,5:7);
EM2_rot_mat = quat2rotm(StretchedSignal{1,2}(:,9:12));
EM2_trans_vec = StretchedSignal{1,2}(:,13:15);

for em_id = 1 : Ref_Length-1
    StretchedSignal{2,1}(em_id+1,:) = FilteredEMGData(em_id+1,:) - FilteredEMGData(em_id,:);
    StretchedSignal{2,2}(em_id+1, 1:6) = [(EM1_trans_vec(em_id+1,:) - EM1_trans_vec(em_id,:)), (EM2_trans_vec(em_id+1,:) - EM2_trans_vec(em_id,:))];
    StretchedSignal{2,2}(em_id+1, 7:8) = [acos((trace(EM1_rot_mat(:,:,em_id) * EM1_rot_mat(:,:,em_id+1)')-1)/2)...
                            acos((trace(EM2_rot_mat(:,:,em_id) * EM2_rot_mat(:,:,em_id+1)')-1)/2)];%*180/pi;
end

StretchedSignal{2, 3}(:, 1) = distMoved(StretchedSignal{1,3}(:, 3))';                                 % ThumbIndexSensor
StretchedSignal{2, 3}(:, 2) = distMoved(StretchedSignal{1,3}(:, 6))';                                 % IndexMiddleSensor 
StretchedSignal{2, 3}(:, 3) = distMoved(StretchedSignal{1,3}(:, 9))';                                 % MiddleRingSensor 
StretchedSignal{2, 3}(:, 4) = distMoved(StretchedSignal{1,3}(:, 12)');                                % RingLittleSensor

StretchedSignal{2, 3}(:, 5) = distMoved(StretchedSignal{1,3}(:, 2),  StretchedSignal{1,3}(:, 5))';    % ThumbIndexEucld
StretchedSignal{2, 3}(:, 6) = distMoved(StretchedSignal{1,3}(:, 5),  StretchedSignal{1,3}(:, 8))';    % IndexMiddleEucld
StretchedSignal{2, 3}(:, 7) = distMoved(StretchedSignal{1,3}(:, 8),  StretchedSignal{1,3}(:, 11))';   % MiddleRingEucld
StretchedSignal{2, 3}(:, 8) = distMoved(StretchedSignal{1,3}(:, 11), StretchedSignal{1,3}(:, 14))';   % RingLittleEucld

StretchedSignal{2, 3}(:, 9) = distMoved(StretchedSignal{1,3}(:, 1),  StretchedSignal{1,3}(:, 2))';    % ThumbNearFarEucld
StretchedSignal{2, 3}(:, 10) = distMoved(StretchedSignal{1,3}(:, 4),  StretchedSignal{1,3}(:, 5))';    % IndexNearFarEucld
StretchedSignal{2, 3}(:, 11) = distMoved(StretchedSignal{1,3}(:, 7),  StretchedSignal{1,3}(:, 8))';    % MiddleNearFarEucld
StretchedSignal{2, 3}(:, 12) = distMoved(StretchedSignal{1,3}(:, 10), StretchedSignal{1,3}(:, 11))';   % RingNearFarEucld
StretchedSignal{2, 3}(:, 13) = distMoved(StretchedSignal{1,3}(:, 13), StretchedSignal{1,3}(:, 14))';   % LittleNearFarEucld

StretchedSignal{2, 4}(:, 1) = double( StretchedSignal{1,4} * 0.0098);
   
%% 4) Generate the data for NN model development
%Fixed-length signal windowing
    WinLenght = 100;  %Every 100 ms
    WinSegments = floor(Ref_Length/WinLenght);    
    n = floor(1+log2(WinLenght));

    %((Ref_Length*8)+(Ref_Length*16)+(Ref_Length*14)+(Ref_Length*1)) %Used to check if WinLenght = 100 is apt
        
    for Wn = 1 : WinSegments
        window = (WinLenght*Wn-WinLenght)+1 : WinLenght*Wn;
        ind = [0 0];
        all_Feats = [];
        for dataId = 1 : length(StretchedSignal)
            ind(1) = ind(2)+1;
            ind(2) = ind(2) + width(StretchedSignal{2,dataId});
            Chunk_To_Process = StretchedSignal{2,dataId}(window(1):window(end), :);
            SpectSegMusc{Cnt+Wn,1}(1:WinLenght, ind(1):ind(2)) = Chunk_To_Process;

            if dataId == 1
                feat_EMG_RMS = rms(Chunk_To_Process);
                feat_ARV = mean(abs(Chunk_To_Process));                  
                all_Feats = [all_Feats, feat_EMG_RMS, feat_ARV];
            else
                all_Feats = [all_Feats, rms(Chunk_To_Process)];                            
            end
        end
        SpectSegMusc{Cnt+Wn,2} = all_Feats;

        if action == "PULL"
            SpectSegMusc{Cnt+Wn,3} = 1;
        elseif action == "PUSH"
            SpectSegMusc{Cnt+Wn,3} = 2;
        end
    end
    
    Cnt = Cnt + Wn;

%Variable-length signals windowing using event-based strategy
%%
figure(44); clf;
LW=2; FsZ=20; axFsZ=18;

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.1], [0.1 0.1], [0.05 0.05]);
%subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

midEMG = floor(size(emg.Data,1)/2);
h=subplot(3,4,1); plot(emg.Data(midEMG-(emg.Resolution/2) : (midEMG+emg.Resolution/2),1),'-r', 'LineWidth',LW);
ylabel('Amplitude', 'FontWeight','bold', 'FontSize',FsZ); 
title({'\fontsize{20}Muscle Activity Signal', ' '});
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'A1', 'Location','southeast'); legend('boxoff');

midEM = floor(size(em.Data,1)/2);
h=subplot(3,4,2); plot(str2double(em.Data(midEM-(em.Resolution/2) : (midEM+em.Resolution/2),1)),'-b', 'LineWidth',LW);
title({'\fontsize{20}Hand Motion Signal', ' '}); ytickformat('%.2f') 
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'A2', 'Location','southeast'); legend('boxoff');

midGLV= floor(size(glv.Data,1)/2);
h=subplot(3,4,3); plot(glv.Data(midGLV-(glv.Resolution/2) : (midGLV+glv.Resolution/2),1),'-k', 'LineWidth',LW);
title({'\fontsize{20}Finger Motion Signal', ' '});  ytickformat('%.3f') 
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'A3', 'Location','southeast'); legend('boxoff');

midFRC = floor(size(force.Data,1)/2);
h=subplot(3,4,4); plot(force.Data(midFRC-(force.Resolution/2) : (midFRC+force.Resolution/2),1),'-g', 'LineWidth',LW);
title({'\fontsize{20}Tool-Vessel Force Data', ' '});   
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'A4', 'Location','southeast'); legend('boxoff');


samp = midEMG-(Ref_resolution/2) : midEMG+(Ref_resolution/2);
h=subplot(3,4,5); plot(sampled_EMG(samp, 1),'-r', 'LineWidth',LW);
ylabel('Amplitude', 'FontWeight','bold', 'FontSize',FsZ); 
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'B1', 'Location','southeast'); legend('boxoff');

h=subplot(3,4,6); plot(sampled_EM(samp, 1),'-b', 'LineWidth',LW);
ytickformat('%.3f');
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'B2', 'Location','southeast'); legend('boxoff');

h=subplot(3,4,7); plot(sampled_glv(samp, 1),'-k', 'LineWidth',LW);
ytickformat('%.3f');
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'B3', 'Location','southeast'); legend('boxoff');

h=subplot(3,4,8); plot(sampled_force(samp, 1),'-g', 'LineWidth',LW);
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'B4', 'Location','southeast'); legend('boxoff');


% samp2 = floor(size(StretchedSignal{1,1},1)/2)-(Ref_resolution/2) : floor(size(StretchedSignal{1,1},1)/2)+(Ref_resolution/2);
% h=subplot(3,4,9); plot(StretchedSignal{1,1}(samp2, 1),'-r', 'LineWidth',LW);
% ylabel('Amplitude', 'FontWeight','bold', 'FontSize',FsZ);
% ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
% 
% h=subplot(3,4,10); plot(StretchedSignal{1,2}(samp2, 1),'-b', 'LineWidth',LW);
% ytickformat('%.3f');
% ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
% 
% h=subplot(3,4,11); plot(StretchedSignal{1,3}(samp2, 1),'-k', 'LineWidth',LW);
% ytickformat('%.3f');
% ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
% 
% h=subplot(3,4,12); plot(StretchedSignal{1,4}(samp2, 1),'-g', 'LineWidth',LW);
% ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';


samp3 = floor(size(StretchedSignal{2,1},1)/2)-(Ref_resolution/2) : floor(size(StretchedSignal{2,1},1)/2)+(Ref_resolution/2);
h=subplot(3,4,9); 
plot(StretchedSignal{2,1}(samp3, 1),'-r', 'LineWidth',LW); 
ylabel('Amplitude', 'FontWeight','bold', 'FontSize',FsZ); 
xlabel('Timestamp (ms)', 'FontWeight','bold', 'FontSize', FsZ);
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'C1', 'Position',[0.152,0.26,0.047,0.03]); legend('boxoff');

h=subplot(3,4,10); 
rescale = MinMax2(StretchedSignal{2,1}(samp3, 1), min(min(StretchedSignal{2,2}(samp3, 1))), max(max(StretchedSignal{2,2}(samp3, 1))));
plot(rescale,'-r', 'LineWidth',LW); hold on;
plot(StretchedSignal{2,2}(samp3, 1),'-b', 'LineWidth',LW); hold off; 
xlabel('Timestamp (ms)', 'FontWeight','bold', 'FontSize', FsZ);
ytickformat('%.3f');
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'C2', 'Position',[0.4,0.26,0.047,0.03]); legend('boxoff');

h=subplot(3,4,11); 
rescale = MinMax2(StretchedSignal{2,1}(samp3, 1), min(min(StretchedSignal{2,3}(samp3, 1))), max(max(StretchedSignal{2,3}(samp3, 1))));
plot(rescale,'-r', 'LineWidth',LW); hold on;
plot(StretchedSignal{2,3}(samp3, 1),'-k', 'LineWidth',LW); hold off; 
xlabel('Timestamp (ms)', 'FontWeight','bold', 'FontSize', FsZ);
ytickformat('%.3f');
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'C3', 'Position',[0.65,0.26,0.047,0.03]); legend('boxoff');

h=subplot(3,4,12); 
rescale = MinMax2(StretchedSignal{2,1}(samp3, 1), min(min(StretchedSignal{2,4}(samp3, 1))), max(max(StretchedSignal{2,4}(samp3, 1))));
plot(rescale,'-r', 'LineWidth',LW); hold on;
plot(StretchedSignal{2,4}(samp3, 1),'-g', 'LineWidth',LW); hold off; 
xlabel('Timestamp (ms)', 'FontWeight','bold', 'FontSize', FsZ);
ax = ancestor(h, 'axes'); ax.FontSize = axFsZ; ax.FontWeight='bold';
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh, 'C4', 'Position',[0.90,0.26,0.047,0.03]); legend('boxoff');


set(gcf,'color','w');