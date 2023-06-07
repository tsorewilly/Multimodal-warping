function FolderFiles = arrangeFiles(directory)
    Folder = dir(directory);
    Folder = Folder(~[Folder.isdir]);
    [~,Idx] = sort([Folder.date]);
    for i = 1 : length(Idx)
        FolderFiles{i,1} = Folder(Idx).name;
    end
end