clear; clc;
warning on;
FolderDir = 'G:\Documents\My Postdoc Work\Intelligent Robot Navigation\EMG\Data\'; 
Lists=dir(fullfile(FolderDir, '*.mat'));
recordLists = {Lists.name};
dataNum= length(recordLists);
st=0;

for redId = 1 : dataNum
    %Load single user's data features
    load([FolderDir, char(recordLists(redId))]);
    fileNameXters = split(char(recordLists(redId)), '-');
    TrialId = split(char(fileNameXters(4)), '.');
    switch(lower(char(fileNameXters{1})))
        case 'omisore'
            UserId = 1;
        case 'duanwenke'
            UserId = 2;
        case 'zhuyisheng'
            UserId = 3;
        case 'wanglei'
            UserId = 4;
        case 'liyifa'
            UserId = 5;
        case 'zhuwang'
            UserId = 6;
        case 'xiongjing'
            UserId = 7;
        otherwise
            UserId = 0;
            warning('Unknown User');
    end
    
    recNum = size(featEMGRMS, 1);
    %T_Id = ones(recNum,1)*str2double(TrialId{1});
    U_Id = ones(recNum,1)*UserId;
    R_Id = ones(recNum,1)*redId;
    
    UserData(st+1:st+recNum, :) = [R_Id, U_Id, featEMGRMS, featEMGMPF, featEMGSamEn, double(featClampForce), featEMinterHand3D,... 
        featEMintraIndex3D, featEMintraThumb3D, featGLinterTIS, featGLinterIMS, featGLinterMRS, featGLinterRLS, featGLinterTIF,...
        featGLinterIMF, featGLinterMRF, featGLinterRLF, featGLintraTNF, featGLintraINF, featGLintraMNF, featGLintraRNF, featGLintraLNF];
    
    st = st+recNum;
end