clear; clc;
%Abductor pollicis brevis (APB), Flexor carpi radialis (FCR), Dorsal interossei (DI), and Extensor carpi radialis (ECR)
frcFolder = 'E:\Rabbit\2020.05.22 (Complete Experiment)\DuanWenke-02\force data\';
emgFolder = 'E:\Rabbit\2020.05.22 (Complete Experiment)\DuanWenke-02\emg data\';
emFolder =  'E:\Rabbit\2020.05.22 (Complete Experiment)\DuanWenke-02\EM data\';
glvFolder = 'E:\Rabbit\2020.05.22 (Complete Experiment)\DuanWenke-02\glove data\';

%% Load Timeseries signals
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

%Load EM Force Data
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

%Load Raw Glove Signal
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

%% Update and Align the signals
%Update timing in the multiple time series signals
dataOut = updateData(em, emg, force, glv);
%Align the updated signals
[em, emg, force, glove] = AlignMultimodalData(dataOut);

vars2keep = {'em', 'emg', 'force', 'glove', 'fileEMG'};
clearvars('-except',vars2keep{:});
%% Process the Timeseries Signal further for Individual Feature Extraction 
figure(2); clf;figure(3); clf;figure(4); clf;figure(5); clf;
figure(1); clf;
    subplot(421);plot(emg.Data(emg.AlignStart : emg.AlignEnd,5)); 
    xlabel('time (msec)');ylabel('EMG value (mV)');title('Abductor pollicis brevis EMG');
    subplot(423);plot(emg.Data(emg.AlignStart : emg.AlignEnd,6)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Flexor carpi radialis EMG');
    subplot(425);plot(emg.Data(emg.AlignStart : emg.AlignEnd,7)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Dorsal interossei EMG');
    subplot(427);plot(emg.Data(emg.AlignStart : emg.AlignEnd,8)); xlabel('time (msec)');ylabel('EMG value (mV)');title('Extensor carpi radialis EMG'); 

    [halfFFT1, offset_signal1, medfreq1] = process_emg(emg.Data(:,1), emg.Resolution); subplot(422);plot(halfFFT1(:,1), halfFFT1(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT2, offset_signal2, medfreq2] = process_emg(emg.Data(:,2), emg.Resolution); subplot(424);plot(halfFFT2(:,1), halfFFT2(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT3, offset_signal3, medfreq3] = process_emg(emg.Data(:,3), emg.Resolution); subplot(426);plot(halfFFT3(:,1), halfFFT3(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 
    [halfFFT4, offset_signal4, medfreq4] = process_emg(emg.Data(:,4), emg.Resolution); subplot(428);plot(halfFFT4(:,1), halfFFT4(:,2)); xlabel('Frequency (HZ)'); ylabel('Signal power (DB)'); 

figure(2);
    [FilteredEMGData(:,1), EMG_Segments(:,:,1), RMS(1,:), Xscale(1,:), rmsXscale(1,:)] = filtEMGSignal(emg.Data(:, 1), emg.Resolution, 50, 500, 64);    
        subplot(411); plot(Xscale(1,:),abs(emg.Data(:, 1)), '-b'); hold on; plot (rmsXscale(1,:), RMS(1,:), '-r', 'linewidth', 2); axis tight
        xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
    [FilteredEMGData(:,2), EMG_Segments(:,:,2), RMS(2,:), Xscale(2,:), rmsXscale(2,:)] = filtEMGSignal(emg.Data(:, 2), emg.Resolution, 50, 500, 64);    
        subplot(412); plot(Xscale(2,:),abs(emg.Data(:, 2)), '-b'); hold on; plot (rmsXscale(2,:), RMS(2,:), '-r', 'linewidth', 2); axis tight
        xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
    [FilteredEMGData(:,3), EMG_Segments(:,:,3), RMS(3,:), Xscale(3,:), rmsXscale(3,:)] = filtEMGSignal(emg.Data(:, 3), emg.Resolution, 50, 500, 64);    
        subplot(413); plot(Xscale(3,:),abs(emg.Data(:, 3)), '-b'); hold on; plot (rmsXscale(3,:), RMS(3,:), '-r', 'linewidth', 2); axis tight
        xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
    [FilteredEMGData(:,4), EMG_Segments(:,:,4), RMS(4,:), Xscale(4,:), rmsXscale(4,:)] = filtEMGSignal(emg.Data(:, 4), emg.Resolution, 50, 500, 64);    
        subplot(414); plot(Xscale(4,:),abs(emg.Data(:, 4)), '-b'); hold on; plot (rmsXscale(4,:), RMS(4,:), '-r', 'linewidth', 2); axis tight
        xlabel('Time in seconds'); ylabel('Amplitude in Millivolts'); title ('Raw Data and its RMS Fit'); legend('Raw Data', 'RMS Fit', 'NumColumns', 2);
%%
%Windowing the EMG signals; 
%An approach with fixed window length and variable window number is adopted
    EMGWinLent = 100;  %Every 2 seconds
    rEMG = floor((emg.AlignEnd-emg.AlignStart)/EMGWinLent);
    n = floor(1+log2(EMGWinLent));
    emg.Data=FilteredEMGData;
    %================================================================feature RMS AVR, MPF, SamEn===============================%
    for i = 1:size(emg.Data,2)
        for j = 1:rEMG
            %disp([a*j-a+1 a*j])
            window = EMGWinLent*j-EMGWinLent+1 : EMGWinLent*j;
            EMGSegId(j,:) = [window(1) window(end)];
            featEMGRMS(j,i) = rms(emg.Data(window, i));
            %featARV(j,i) = mean(abs(emg.Data(window, i)));  
            featEMGMPF(j,i) = mpf1(emg.Data(window, i),1000);
            Std(j,i) = std(emg.Data(window, i));  
            AllSamEn = sampenc1(emg.Data(window, i),2,0.2*Std(j,i));
            featEMGSamEn(j,i) = AllSamEn(1);            

            [counts, centers] = hist(emg.Data(window, i),n);        % Calculate the histograms            
            hdat(j,:) = counts./sum(counts);                        % Convert histograms to probability values            
            featEMGEntroy(j,i) = -sum(hdat(j,:).*log2(hdat(j,:)+eps)); % Compute entroy
        end
    end

%%
%Window the Force signals using equivalent of EMG Segmentation indices; 
    ForceSegId = floor(MinMax(EMGSegId, max(max(EMGSegId)), min(min(EMGSegId)), force.AlignStart, force.AlignEnd));
    for i = 1:rEMG%floor(window(end)/size(fVals,1))%Extract segments in same time domain in EMG
        featClampForce(i,1) = double(rms(force.Data(ForceSegId(i,1):ForceSegId(i,2)))*0.0098);
    end

%%
    ThumbFinger = deNoiseEM(double(em.Data(:, 1:7)));
    IndexFinger = deNoiseEM(double(em.Data(:, 9:15)));
    interHand3D = distMoved(ThumbFinger(:,5:7), IndexFinger(:,5:7))';
    intraThumb3D = distMoved(ThumbFinger(:,5:7))';
    intraIndex3D = distMoved(IndexFinger(:,5:7))';

%Window the EM signals using equivalent of EMG Segmentation indices; 
    EMSegId = floor(MinMax(EMGSegId, max(max(EMGSegId)), min(min(EMGSegId)), em.AlignStart, em.AlignEnd));
    for i = 1:rEMG%floor(window(end)/size(emData,1)-1)-1%Extract segments in same time domain in EMG
        featEMinterHand3D(i,1) = rms(interHand3D(EMSegId(i, 1):EMSegId(i, 2)));
        featEMintraThumb3D(i,1) = rms(intraThumb3D(EMSegId(i, 1):EMSegId(i, 2)));
        featEMintraIndex3D(i,1) = rms(intraIndex3D(EMSegId(i, 1):EMSegId(i, 2)));
    end
%%
%Window the GLove signals using equivalent of EMG Segmentation indices; 
%This features are planar unlike the EM features
    GLSegId = floor(MinMax(EMGSegId, max(max(EMGSegId)), min(min(EMGSegId)), glove.AlignStart, glove.AlignEnd));
    ThumbIndexSensor = distMoved(glove.Data(:, 3))'; 
    IndexMiddleSensor = distMoved(glove.Data(:, 6)); 
    MiddleRingSensor = distMoved(glove.Data(:, 9)); 
    RingLittleSensor = distMoved(glove.Data(:, 12)); 

    ThumbIndexEucld = distMoved(glove.Data(:, 2), glove.Data(:, 5))'; 
    IndexMiddleEucld = distMoved(glove.Data(:, 5), glove.Data(:, 8))'; 
    MiddleRingEucld = distMoved(glove.Data(:, 8), glove.Data(:, 11))'; 
    RingLittleEucld = distMoved(glove.Data(:, 11), glove.Data(:, 14))'; 

    ThumbNearFarEucld = distMoved(glove.Data(:, 1), glove.Data(:, 2))'; 
    IndexNearFarEucld = distMoved(glove.Data(:, 4), glove.Data(:, 5))'; 
    MiddleNearFarEucld = distMoved(glove.Data(:, 7), glove.Data(:, 8))'; 
    RingNearFarEucld = distMoved(glove.Data(:, 10), glove.Data(:, 11))'; 
    LittleNearFarEucld = distMoved(glove.Data(:, 13), glove.Data(:, 14))'; 

        for i = 1:rEMG
        %Compute features based interfinger sensors
        featGLinterTIS(i,1) = rms(ThumbIndexSensor(GLSegId(i, 1):GLSegId(i, 2)));       %TIS: ThumbIndex Sensor
        featGLinterIMS(i,1) = rms(IndexMiddleSensor(GLSegId(i, 1):GLSegId(i, 2)));      %IMS: IndexMiddle Sensor
        featGLinterMRS(i,1) = rms(MiddleRingSensor(GLSegId(i, 1):GLSegId(i, 2)));       %TRS: MiddleRing Sensor
        featGLinterRLS(i,1) = rms(RingLittleSensor(GLSegId(i, 1):GLSegId(i, 2)));       %RLS: RingLittle Sensor

        %Add more features from interfinger Euclidean movements
        featGLinterTIF(i,1) = rms(ThumbIndexEucld(GLSegId(i, 1):GLSegId(i, 2)));       %TIF: ThumbIndex Far
        featGLinterIMF(i,1) = rms(IndexMiddleEucld(GLSegId(i, 1):GLSegId(i, 2)));       %TMF: IndexMiddle Far
        featGLinterMRF(i,1) = rms(MiddleRingEucld(GLSegId(i, 1):GLSegId(i, 2)));       %TRF: MiddleRing Far
        featGLinterRLF(i,1) = rms(RingLittleEucld(GLSegId(i, 1):GLSegId(i, 2)));       %IMF: RingLittle Far

        %Additional features from intrafinger Euclidean movements 
        %Compute at two consecutive points
        featGLintraTNF(i,1) = rms(ThumbNearFarEucld(GLSegId(i, 1):GLSegId(i, 2)));     %TNF: Thumb NearFar
        featGLintraINF(i,1) = rms(IndexNearFarEucld(GLSegId(i, 1):GLSegId(i, 2)));     %INF: Index NearFar
        featGLintraMNF(i,1) = rms(MiddleNearFarEucld(GLSegId(i, 1):GLSegId(i, 2)));     %MNF: Middle NearFar
        featGLintraRNF(i,1) = rms(RingNearFarEucld(GLSegId(i, 1):GLSegId(i, 2)));     %RNF: Ring NearFar
        featGLintraLNF(i,1) = rms(LittleNearFarEucld(GLSegId(i, 1):GLSegId(i, 2)));     %LNF: Little NearFar
    end

%%
%save data
cathFile = extractBefore(deblank(fileEMG), length(deblank(fileEMG))-3);    
save (['Data/',cathFile,'.mat'], 'featEMGRMS', 'featEMGMPF', 'featEMGSamEn', 'featClampForce',...
    'featEMinterHand3D', 'featEMintraIndex3D', 'featEMintraThumb3D', 'featGLinterTIS', 'featGLinterIMS', ...
    'featGLinterMRS', 'featGLinterRLS', 'featGLinterTIF', 'featGLinterIMF', 'featGLinterMRF', 'featGLinterRLF',...
    'featGLintraTNF', 'featGLintraINF', 'featGLintraMNF', 'featGLintraRNF', 'featGLintraLNF');

%%
%Display processed sigals
%if(fileEMG==0 || fileForce==0)
    figure(3);clf
    subplot(411); plot(featEMGRMS(:,1),'-r*', 'LineWidth',3); hold on; plot(featEMGRMS(:,2),'-b*', 'LineWidth',3);  plot(featEMGRMS(:,3),'-k*', 'LineWidth',3);  plot(featEMGRMS(:,4),'-g*', 'LineWidth',3); hold off;
    ylabel('Root mean square');     xlabel('Window Index'); title('RMS Analysis Window');   legend('APB','FCR','DI','ECR', 'Location','best', 'NumColumns', 2);
    subplot(412); plot(featEMGSamEn(:,1),'-r*', 'LineWidth',3); hold on; plot(featEMGSamEn(:,2),'-b*', 'LineWidth',3);  plot(featEMGSamEn(:,3),'-k*', 'LineWidth',3);  plot(featEMGSamEn(:,4),'-g*', 'LineWidth',3); hold off;
    ylabel('Sample Entropy');     xlabel('Window Index');   title('Sample Entroy Analysis Window');  legend('APB','FCR','DI','ECR', 'Location','best', 'NumColumns', 2);
%     subplot(412); plot(featARV(:,1),'-r*'); hold on; plot(featARV(:,2),'-b*');  plot(featARV(:,3),'-k*');  plot(featARV(:,4),'-g*'); hold off;
%     ylabel('Sample Entropy');     xlabel('Window Index');   title('Sample Entroy Analysis Window');  legend('APB','FCR','DI','ECR', 'Location','best', 'NumColumns', 2);
    subplot(413); plot(featEMGMPF(:,1),'-r*', 'LineWidth',3); hold on; plot(featEMGMPF(:,2),'-b*', 'LineWidth',3);  plot(featEMGMPF(:,3),'-k*', 'LineWidth',3);  plot(featEMGMPF(:,4),'-g*', 'LineWidth',3); hold off;
    ylabel('Mean Power Frequency (Hz)');     xlabel('Window Index'); title('MPF Analysis Window');  legend('APB','FCR','DI','ECR', 'Location','best', 'NumColumns', 2);
    subplot(414); plot(featEMGEntroy(:,1),'-r*', 'LineWidth',3); hold on; plot(featEMGEntroy(:,2),'-b*', 'LineWidth',3);  plot(featEMGEntroy(:,3),'-k*', 'LineWidth',3);  plot(featEMGEntroy(:,4),'-g*', 'LineWidth',3); hold off;
    ylabel('Feature Entroy');     xlabel('Bin Index');      title('Entroy Analysis Window');    legend('APB','FCR','DI','ECR', 'Location','best', 'NumColumns', 2);

    figure(4);clf
    subplot(2,2,[1 2]); plot(force.Data, '-r', 'LineWidth',3); ylabel('Proximal Force (gram)');     xlabel('Time'); title('Analysis of Catheterization Force');
    subplot(223); plot(featClampForce,'-r*', 'LineWidth',3);
    subplot(224); histfit(force.Data)

    figure(5); clf    
    plot(featEMinterHand3D, '-ro', 'LineWidth',3); hold on;     plot(featEMintraIndex3D,'-k*', 'LineWidth',3);    plot(featEMintraThumb3D,'-bp', 'LineWidth',3);
%    subplot(224); histfit(fVals)
    ylabel('Spatial Translational Motion(mm)');     xlabel('Time');    title('Hand and Finger Movements during Chateterization');   
     legend('Hand Movement','Thumb Movement','Index Movement', 'Location','best');
    hold off;

    figure(6);clf
    subplot(311); plot(featGLinterTIS,'-r*', 'LineWidth',3); hold on; plot(featGLinterIMS,'-b*', 'LineWidth',3);  plot(featGLinterMRS,'-k*', 'LineWidth',3);  plot(featGLinterRLS,'-g*', 'LineWidth',3); hold off;
    ylabel('RMS Values');     xlabel('Window Index'); title('Inter-Finger Sensor Data');   legend('TIS','TMS','TRS','RLS', 'Location','best', 'NumColumns', 2);
    subplot(312); plot(featGLinterTIF,'-r*', 'LineWidth',3); hold on; plot(featGLinterIMF,'-b*', 'LineWidth',3);  plot(featGLinterMRF,'-k*', 'LineWidth',3);  plot(featGLinterRLF,'-g*', 'LineWidth',3); hold off;
    ylabel('RMS Values');     xlabel('Window Index'); title('Inter-Finger Euclidean Distance');   legend('TIF','IMF','TRF','IMF', 'Location','best', 'NumColumns', 2);
    subplot(313); plot(featGLintraTNF,'-r*', 'LineWidth',3); hold on; plot(featGLintraINF,'-b*', 'LineWidth',3);  plot(featGLintraMNF,'-k*', 'LineWidth',3);  plot(featGLintraRNF,'-g*', 'LineWidth',3); plot(featGLintraLNF,'-c*', 'LineWidth',3); hold off;
    ylabel('RMS Values');     xlabel('Window Index'); title('Intra-Finger Euclidean Distance');  legend('TNF','INF','MNF','RNF','LNF', 'Location','best', 'NumColumns', 2);

%end    
