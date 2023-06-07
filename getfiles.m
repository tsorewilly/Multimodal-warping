function fileInfo = getfiles(filename)
% This function search for the numbre of files written or created on 
% specific date and plot a bargraph. The filename is a string variable as
% below.
% filename = 'S:\Support Operations\Support Data\SRS MapCHECK';
% Then, type in Matlab command prompt as below
% >> getfileinfo('S:\SO\Support Data\Product1')
% This function is not idiot proof - Erros can occur if used outside of
% scope.
tic
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
        dateCreated(i,1) = datetime(cell2mat(strcat(c{1}{k-2},{' '},c{1}{k-1},{' '},c{1}{k})),'InputFormat','MM/dd/uuuu hh:mm aa');
        n = length(str2num(c{1}{k+1}));
        fileSize(i,1) = sum(fliplr(str2num(c{1}{k+1})).*(10.^(3.*((0:n-1)))));
        i = i + 1;
    end
end
% Sort according to date of creation
[dateCreated,id] = sort(dateCreated);
fileSize = fileSize(id);
% plot a bar graph
span = days(dateCreated(end)-dateCreated(1));
t = span/(length(fileSize)-1);
figure('Name',filename,'MenuBar', 'none','ToolBar', 'none')
if ~isequal(span, 0)
    bar(0:t:span,cumsum(fileSize)/2^20,'stacked')
    xlabel('Duratiion (days)')
else
    bar(cumsum(fileSize)/2^20,'stacked')
    xlabel('No. of Files over 1 day')
end
    ylabel('Cumulative data (MB)')
    title(filename)
info = sprintf('Contains: %d files, %d, bytes.\nRate of Data Accum. ~ %0.2G GB/year.\nWritten dates:\n  Start: %s\n  End: %s',...
    length(fileSize), sum(fileSize),(365*sum(fileSize))/(2^30*span),...
    dateCreated(1),dateCreated(end));
text(span*0.02,sum(fileSize)*.9/2^20,info);
toc
end
