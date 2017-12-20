%% 1. PLOT OF STATIONS IN THE PHILIPPINES

clear all; close all;

% read stations and coordinates file
T = readtable('./data_slv/station_coords.txt', 'Delimiter',',')

h1=figure();
% Map of the Philippines
xlims=[116 128];    % longitude limits (deg E)
ylims=[4 20];       % latitude limits (deg N)
m_proj('mercator','lon',xlims,'lat',ylims); % map projection
m_gshhs_i()         % Generate coastline
m_grid
hold on
point_color = {'w','r'}
for i=1:length(T.Longitude)
    m_plot(T.Longitude(i),T.Latitude(i),'ro','MarkerEdgeColor','k',...
        'MarkerFaceColor',point_color{T.Complete(i)+1})
end
m_text(T.Longitude+0.2,T.Latitude,T.Station)

% % Save figure to file dialog box
% display(sprintf('Saving figure to file...'))
% [filename,pathname]=uiputfile('*.png','Save map of statiions as (*.png)');
% saveas(h1,[pathname filename]);
% display(sprintf('Figure saved to %s%s',pathname,filename))

%% 2. COMPILE DATA - HOURLY

clear all; close all;
% read description of stations, datasets, and filenames
T = readtable('./data_slv/filenames_description.txt', 'Delimiter',',');

% hourly
rows = find(strcmp(T.Resolution,'hourly'));
names = T.Varname(rows);
figure()
for i=1:length(rows)
    fname = strcat('./data_slv/',T.Filename{rows(i)});
    tmp = csvread(fname);
    year = tmp(:,1); month = tmp(:,2); day = tmp(:,3); hour = tmp(:,4);
    minute = repmat(0,size(tmp(:,1)));
    second = minute;
    time = datenum(year,month,day,hour,minute,second);
    % set flagged rows to NaN values
    flaggedrows = find(isnan(tmp(:,5)));
    tmp(flaggedrows,5) = NaN;
    hourly(i).name = names{i};
    hourly(i).raw = horzcat(time,tmp(:,5));
    hourly(i).lat = T.Latitude(rows(i));
    hourly(i).lon = T.Longitude(rows(i));
    % plot all time series for inspection
    subplot(6,2,i)
    plot(time,tmp(:,5),'k-')
    xlim([time(1) time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',hourly(i).name,...
        ' (',num2str(hourly(i).lat),...
        ', ',num2str(hourly(i).lon),')'])
end

save('./output/hourly.mat','hourly');

% change to save to output folder
%% 3. PLOT POWER SPECTRA OF ALL STATIONS

clear all; close all;
load('./output/hourly.mat')

for i = 1:length(hourly)
    fig_powerspectra(hourly(i).raw(:,2),hourly(i).raw(:,1),hourly(i).lat);
end

%% 4. TIDAL COMPONENTS ANALYSIS WITH t_tide() 

clear all; close all;
load('./output/hourly.mat')

% create a new series, regularly gridded time, 
% set unsampled times to NaN
for i=1:length(hourly)
    tmp = hourly(i).raw;
    dt = min(diff(tmp(:,1)));
    % regular grid
    tvec = (tmp(1,1):dt:tmp(end,1))';
    nt = length(tvec);
    xvec = NaN(nt,1);
    [m,locb] = ismember(round(tmp(:,1),4),round(tvec,4));
    xvec(locb) = tmp(:,2);
    xmu = nanmean(xvec);
    % regular time grid
    hourly(i).time = tvec;
    % observed slv at regular time grid
    hourly(i).slv = xvec;
    % convert observed slv to anomaly values
    hourly(i).aslv = xvec-repmat(xmu,length(xvec),1);
end

% tide analysis
tidemat = [];
for i=1:length(hourly)
    % need to truncate long time series with many gaps ('mnla', 'leg')
    % so t_tide will work
    maxobs = 163000;
    if(strcmp(hourly(i).name,'mnla')|strcmp(hourly(i).name,'leg'))
        hourly(i).tides = calc_tides(hourly(i).aslv(1:maxobs),...
            hourly(i).time(1:maxobs),hourly(i).lat);
    else
        hourly(i).tides = calc_tides(hourly(i).aslv,...
            hourly(i).time,hourly(i).lat);
    end
    tidemat = vertcat(tidemat,[hourly(i).lon,...
        hourly(i).lat,hourly(i).tides.relcon]);
    % check tide predictions and observed sea level anomalies
    figure()
    plot(hourly(i).time,hourly(i).aslv,'k-'); hold on
    plot(hourly(i).tides.T,hourly(i).tides.tide,'r')
    xlim([hourly(i).time(1) hourly(i).time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    
end

save('./output/hourly_analysis.mat','hourly','tidemat');


%% 4. PLOT SPATIAL VARIABILITY IN THE RELATIVE IMPORTANCE OF DIURNAL AND
% SEMI-DIURNAL COMPONENTS

clear all; close all;
load('./output/hourly_analysis.mat')

h1=figure();
% Map of the Philippines
xlims=[116 128];    % longitude limits (deg E)
ylims=[4 20];       % latitude limits (deg N)
m_proj('mercator','lon',xlims,'lat',ylims); % map projection
m_gshhs_i()         % Generate coastline
m_grid
hold on

% Take only a subset since 2 stations have 2 analyses;
% Retain only the analysis based on most recent data given station
selectrows = [2,3,5:11];
x = tidemat(selectrows,1);
y = tidemat(selectrows,2);
z = tidemat(selectrows,3);

[xi, yi] = meshgrid(...
    linspace(min(x),max(x)),...
    linspace(min(y),max(y)));

zi = griddata(x,y,z, xi,yi);
%m_contourf(xi,yi,zi)
m_pcolor(xi,yi,zi); shading flat; colorbar
h = colorbar; ylabel(h,'(M_2+S_2) / (O_1+K_1)')

% Save figure to file dialog box
display(sprintf('Saving figure to file...'))
[filename,pathname]=uiputfile('*.png','Save map of statiions as (*.png)');
saveas(h1,[pathname filename]);
display(sprintf('Figure saved to %s%s',pathname,filename))


