function fileSize = get2fileinfo(filename)
    % This function search for  of files written or created on specific date and plot a bargraph.      
    % DOS string for "FILE WRITTEN" data indicating /tw. 
    % For 'File Creation Date', use /tc
    dosStr = char(strcat({'dir /s /tw '}, {'"'}, filename,{'"'}));
    % Reading DOS results
    [~,results] = dos(dosStr);
    c = textscan(results,'%s');
    % Extract file size (fileSize) and date created (dateCreated)
    i = 1;
    for k = 1:length(c{1})
        if (isequal(c{1}{k},'PM')||isequal(c{1}{k},'AM'))&&(~isequal(c{1}{k+1},'<DIR>'))
            dateCreated = datetime(cell2mat(strcat(c{1}{k-2},{' '},c{1}{k-1},{' '},c{1}{k})),'Format','yyyy-MM-dd HH:mm:ss.SSS');
        end
    end
end

