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

%% 1)	Normalization and Selection of Reference signal
%Update timing in the multiple time series signals
%dataOut = updateData(em, emg, force, glv); Not used here, see
%'ManualFeatureExtractor.m'

%resample data to same resolution frequency 

Ref_resolution = emg.Resolution; %(Set reference as EM)
sampled_EMG = resample(emg.Data, Ref_resolution, emg.Resolution);
sampled_EM = resample(cellfun(@(x)str2double(x), em.Data), Ref_resolution, em.Resolution);
sampled_EM(sampled_EM<-200) = -200; sampled_EM(sampled_EM>200) = 200; 
sampled_glv = resample(glv.Data, Ref_resolution, glv.Resolution);
sampled_force = resample(double(force.Data), Ref_resolution, force.Resolution);

figure(1)
subplot(4,2,[1,2]);plot(sampled_EMG(:,1)); 
subplot(4,2,[3,4]);plot(sampled_EM(:,1)); 
subplot(4,2,[5,6]);plot(sampled_glv(:,1)); 
subplot(4,2,[7,8]);plot(sampled_force(:,1));

%Window the the reference signal; 
%An approach with fixed window length and variable window number is adopted
n = 20; %Every 2 seconds

[TR_rEMG, VR_rEMG] = windowSignal(sampled_EMG(:,1), n);
[TR_rEM, VR_rEM] = windowSignal(sampled_EM(:,1), n);
[TR_rglv, VR_rglv] = windowSignal(sampled_glv(:,1), n);
[TR_rforce, VR_rforce] = windowSignal(sampled_force(:,1), n);

figure(2)
subplot(4,2,[1,2]);plot(VR_rEMG'); 
subplot(4,2,[3,4]);plot(VR_rEM'); 
subplot(4,2,[5,6]);plot(VR_rglv'); 
subplot(4,2,[7,8]);plot(VR_rforce');

%Signal Referencing and Stretching
[LongestSeq, Id] = max([length(VR_rEMG), length(VR_rEM), length(VR_rglv), length(VR_rforce)]);
WWL = 0.5; %WWL
%Make signals equal length, and fill their contents from both rears
[EL_rEMG, EL_rEM, EL_rglv, EL_rforce] = deal(zeros(n, LongestSeq));
%methods = [1:'inverse_proportion'; 2: 'inverse_squared']
for i = 1:n
    for j = 1:floor(LongestSeq/2)
        EL_rEMG(i,j) = VR_rEMG(i,j); EL_rEMG(i,end-j+1) = VR_rEMG(i,end-j+1);
        EL_rEM(i,j) = VR_rEM(i,j); EL_rEM(i,end-j+1) = VR_rEM(i,end-j+1);
        EL_rglv(i,j) = VR_rglv(i,j); EL_rglv(i,end-j+1) = VR_rglv(i,end-j+1);
        EL_rforce(i,j) = VR_rforce(i,j); EL_rforce(i,end-j+1) = VR_rforce(i,end-j+1);       
    end
    %fill median index in case of odd-lenght signal
    if(mod(LongestSeq,2)==1)
        EL_rEMG(i,ceil(LongestSeq/2)) = mean(VR_rEMG(i,j));
        EL_rEM(i,ceil(LongestSeq/2)) = mean(VR_rEM(i,j));
        EL_rglv(i,ceil(LongestSeq/2)) = mean(VR_rglv(i,j));    
        EL_rforce(i,ceil(LongestSeq/2)) = mean(VR_rforce(i,j));
    end
    
    ELD_rEMG(i,:) = EL_rEMG(i,:) - DistShiftSeq(EL_rEMG(i,:), WWL, 1)';
    ELD_rEM(i,:) = EL_rEM(i,:) - DistShiftSeq(EL_rEM(i,:), WWL, 1)'; 
    ELD_rglv(i,:) = EL_rglv(i,:) - DistShiftSeq(EL_rglv(i,:), WWL, 1)'; 
    ELD_rforce(i,:) = EL_rforce(i,:) - DistShiftSeq(EL_rforce(i,:), WWL, 1)'; 
   
end


figure(3)
subplot(4,2,[1,2]);plot(EL_rEMG');     subplot(4,2,[3,4]);plot(EL_rEM'); 
subplot(4,2,[5,6]);plot(EL_rglv');     subplot(4,2,[7,8]);plot(EL_rforce');   

figure(4)
subplot(4,2,[1,2]);plot(ELD_rEMG');     subplot(4,2,[3,4]);plot(ELD_rEM'); 
subplot(4,2,[5,6]);plot(ELD_rglv');     subplot(4,2,[7,8]);plot(ELD_rforce');


%Generating the Warped Signal



%% B2)	Companion reference signal generation

%Align the updated signals
%[em, emg, force, glove] = staticTimeWarping([em, emg, force, glove]); %Warping and Alingment of Multimodal Data done together
%[em, emg, force, glove] = DyanmicTimeWarping([em, emg, force, glove]); %with resolution-based DTW method
dataIn = [em, emg, force, glv];

