function data = deNoiseEM(ddata, amp)
    if(~exist('w','var'))
        amp = 1000; %Use default warping window
    end
    data=ddata;
    data(abs(data)>amp) = NaN;
    for i = 1: size(data,2)
        cleanData = data(:,i);
        cleanData(isnan(cleanData)) = mean(cleanData,'omitnan');
        data(:,i) = cleanData;
    end
%     scatter3(data(:,5),data(:,6),data(:,7),'ro'); hold on;
%     plot3(data(1,5),data(1,6),data(1,7),'bo'); hold off;
end