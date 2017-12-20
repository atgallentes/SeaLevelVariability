%% 1. COMPILE DATA - ANNUAL

clear all; close all;
% read description of stations, datasets, and filenames
T = readtable('/data_slv/filenames_description.txt', 'Delimiter',',');

rows = find(strcmp(T.Resolution,'annual'));
names = T.Varname(rows);

figure()
for i=1:length(rows)
    fname = strcat('/data_slv/',T.Filename{rows(i)});
    tmp = readtable(fname,'Delimiter',';','ReadVariableNames',false); 
    annual(i).name = names{i};
    annual(i).raw = tmp;
    annual(i).lat = T.Latitude(rows(i));
    annual(i).lon = T.Longitude(rows(i));
    % convert flagged rows to NaN;
    flaggedrows = find(tmp.Var2==-99999);
    tmp.Var2(flaggedrows) = NaN;
    annual(i).time = tmp.Var1;
    annual(i).slv = tmp.Var2;
    annual(i).aslv = tmp.Var2 - repmat(nanmean(tmp.Var2),...
        length(tmp.Var2),1);
    
    subplot(6,2,i)
    plot(annual(i).time,annual(i).aslv,'k-')
    xlim([annual(i).time(1) annual(i).time(end)])
    %datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',annual(i).name,...
        ' (',num2str(annual(i).lat),...
        ', ',num2str(annual(i).lon),')'])
end

save('./output/annual.mat','annual');

% To do: mnlb - annual?, 
% Remove from annual: zam, cag
% Calculate for cur, lub? 
% See AVISO data