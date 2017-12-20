%% 1. COMPILE DATA - DAILY

clear all; close all;
% read description of stations, datasets, and filenames
T = readtable('/data_slv/filenames_description.txt', 'Delimiter',',');

rows = find(strcmp(T.Resolution,'daily'));
names = T.Varname(rows);
for i=1:length(rows)
    fname = strcat('/data_slv/',T.Filename{rows(i)});
    tmp = csvread(fname); 
    time = datenum(tmp(:,1:3));   
    daily(i).name = names{i};
    daily(i).raw = horzcat(time,tmp(:,4));
    daily(i).lat = T.Latitude(rows(i));
    daily(i).lon = T.Longitude(rows(i));
    % plot all time series for inspection
    subplot(6,2,i)
    plot(time,tmp(:,4),'k-')
    xlim([time(1) time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',daily(i).name,...
        ' (',num2str(daily(i).lat),...
        ', ',num2str(daily(i).lon),')'])
end



