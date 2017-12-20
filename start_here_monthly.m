%% 1. COMPILE MONTHLY DATA

clear all; close all;
% read description of stations, datasets, and filenames
T = readtable('./data_slv/filenames_description.txt', 'Delimiter',',');

rows = find(strcmp(T.Resolution,'monthly'));
names = T.Varname(rows);
for i=1:length(rows)
    fname = strcat('./data_slv/',T.Filename{rows(i)});    
    tmp = dlmread(fname,';'); 
    %convert year-month decimal to datenum
    %year-month-decimal-form = year + (month-0.5)/12.0  
    year = floor(tmp(:,1));
    decpart = tmp(:,1) - year;
    month = round(decpart*12 + 0.5);
    time = datenum(year,month,repmat(1,size(month)));
    tmp = horzcat(time,tmp(:,2:end));
    monthly(i).name = names{i};
    monthly(i).raw = tmp;
    monthly(i).lat = T.Latitude(rows(i));
    monthly(i).lon = T.Longitude(rows(i));
    % plot all time series for inspection
    subplot(9,2,i)
    plot(time,tmp(:,2),'k-')
    xlim([time(1) time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',monthly(i).name,...
        ' (',num2str(monthly(i).lat),...
        ', ',num2str(monthly(i).lon),')'])
end

% Find flagged values = -99999, set to NaN

figure()
for i=1:length(rows)
    tmp = monthly(i).raw;
    flaggedrows = find(tmp(:,2)==-99999);
    tmp(flaggedrows,2) = NaN;   
    xmu = nanmean(tmp(:,2));
    monthly(i).time = tmp(:,1);
    monthly(i).slv = tmp(:,2);
    % convert observed slv to anomaly values
    monthly(i).aslv = tmp(:,2)-repmat(xmu,length(tmp(:,2)),1);
    
    %check if anomaly values are correct
    subplot(9,2,i)
    plot(monthly(i).time,monthly(i).aslv,'k-')
    xlim([monthly(i).time(1) monthly(i).time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',monthly(i).name,...
        ' (',num2str(monthly(i).lat),...
        ', ',num2str(monthly(i).lon),')'])
end


save('./output/monthly.mat','monthly');

%% 2. PROCESS SOI DATA
clear all; close all;
fname = '/data_params/SOI.txt';

soi = readtable(fname, 'Delimiter',',');
soi.year = floor(soi.Date/100);
soi.month = soi.Date - soi.year*100;
soi.time = datenum(soi.year,soi.month,repmat(1,length(soi.year),1));

save('./output/soi.mat','soi');

%% 3. DETREND DATA - LONG TERM TREND

clear all; close all;
load('./output/monthly.mat')

for i=1:length(monthly)
    statname = monthly(i).name;
    pred_aslv = [];      %initialize array for predicted values
    
    if strcmp(statname,'mnla')   %fit piecewise linear to mnla data
        % trend 
        t1_ind = find(monthly(i).time <= datenum(1966,1,1));
        t2_ind = find(monthly(i).time > datenum(1966,1,1));
        tvec = {t1_ind t2_ind}
            
        for j = 1:2
            lm = regress(monthly(i).aslv(tvec{j}),...
                [repmat(1,length(tvec{j}),1) monthly(i).time(tvec{j})])
            y = repmat(lm(1),length(tvec{j}),1) + lm(2) * monthly(i).time(tvec{j});
            pred_aslv = vertcat(pred_aslv, y) 
            %det_aislv1 = monthly(i).aslv(t1_ind) - y1;
        end                     
    else % fit linear trend
        lm = regress(monthly(i).aslv,...
            [repmat(1,length(monthly(i).aslv),1) monthly(i).time])
        y = repmat(lm(1),length(monthly(i).aslv),1) + lm(2) * monthly(i).time;
        pred_aslv = vertcat(pred_aslv, y) 
    end
    
    % detrend time series and save to variable
    monthly(i).longtrend = pred_aslv;
    monthly(i).aslv_det1 = monthly(i).aslv - pred_aslv;
    
    % check if detrending is correct
    figure()
    plot(monthly(i).time,monthly(i).aslv)
    xlim([monthly(i).time(1) monthly(i).time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',monthly(i).name,...
        ' (',num2str(monthly(i).lat),...
        ', ',num2str(monthly(i).lon),')'])
    hold on;
    plot(monthly(i).time,pred_aslv,'ro',...
        'MarkerSize',.8,'MarkerFaceColor','r');
    plot(monthly(i).time,monthly(i).aslv_det1,'k-',...
        'MarkerSize',.8,'MarkerFaceColor','r');    
end

save('./output/monthly2.mat','monthly')

%% Plot long-term trends

clear all; close all;

load('./output/monthly2.mat')

selectrows = [1,4,5,6,8,9,10,14,17];
tmp = monthly(selectrows)
for i=1:length(tmp)
    subplot(5,2,i)
    plot(tmp(i).time,tmp(i).aslv,'k-'); hold on
    plot(tmp(i).time,tmp(i).longtrend,'ro',...
        'MarkerSize',0.8,'MarkerFaceColor','r');
    xlim([tmp(i).time(1) tmp(i).time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',monthly(i).name,...
        ' (',num2str(monthly(i).lat),...
        ', ',num2str(monthly(i).lon),')'])
    if (mod(i,2) == 1)
        ylabel('aslv (mm)')
    end
    %xlabel('time')
end

%% 4. FURTHER DETREND DATA - ANNUAL AND SEMI-ANNUAL COMPONENTS

clear all; close all;
load('./output/monthly2.mat')

for i=1:length(monthly)
    t = monthly(i).time;
    X = ones(length(t),5);
    w = 2*pi / 12;
    
    for K = 1:length(t)
        wdt = w*K
        wdt2 = wdt*2;
        X(K,1) = 1; %mean
        X(K,2) = cos(wdt);
        X(K,3) = sin(wdt);
        X(K,4) = cos(wdt2);
        X(K,5) = sin(wdt2);
    end
    alpha = 0.05;
    [b,int] = regress(monthly(i).aslv_det1,X,alpha);
    annslv = X(:,1) * b(1) + X(:,2) * b(2) + X(:,3)*b(3) + X(:,4)*b(4) + X(:,5)*b(5);
    monthly(i).aslv_det2 = monthly(i).aslv_det1 - annslv;
    
    % check if detrending is correct
    figure()
    plot(monthly(i).time,monthly(i).aslv_det1)
    xlim([monthly(i).time(1) monthly(i).time(end)])
    datetick('x','mmm-yy','keeplimits','keepticks')
    title(['Station: ',monthly(i).name,...
        ' (',num2str(monthly(i).lat),...
        ', ',num2str(monthly(i).lon),')'])
    hold on;
%     plot(monthly(i).time,annslv,'r--',...
%         'MarkerSize',.8,'MarkerFaceColor','r');
    plot(monthly(i).time,monthly(i).aslv_det2,'k-',...
        'MarkerSize',.8,'MarkerFaceColor','r');
    
end

save('./output/monthly3.mat','monthly')


%% 5. LAG REGRESSION DETRENDED ASLV WITH SOI

clear all; close all;
load('./output/monthly3.mat')   %observations
load('./output/soi.mat')        %SOI

soireg = [];
for i=1:length(monthly)
    t_obs = monthly(i).time;
    ind = find(t_obs >= min(soi.time) & t_obs <= max(soi.time));    
    [m,locb] = ismember(t_obs(ind),soi.time);
    x = monthly(i).aslv_det2(ind);
    y = soi.Value(locb);
    x = naninterp(monthly(i).aslv_det2(ind));
    [z,lag] = conn(x,soi.Value(locb),10); 
    % find out at which lag correlation is highest
    [mval,mind] = max(abs(z));
    optlag = lag(mind);
    
    soireg = vertcat(soireg, [monthly(i).lon,monthly(i).lat,...
        optlag,z(mind)]);
       
% %     figure()
% %     plot(lag,z)
% %     xlabel('lag (months)')
% %     ylabel('correlation')
% %     % note: conn(x,y,n)
% %     % x leads for + lag; y leads for - lag
% %     title([monthly(i).name,' opt lag: ',num2str(optlag),...
% %         '; %cor: ',num2str(mval)]);   
    
% %     figure()
% %     plot(t_obs,monthly(i).aslv_det2,'k-'); hold on
% %     plot(soi.time,soi.Value,'r-')
% %     datetick('x','mmm-yy','keeplimits','keepticks')
% %     title(['Station: ',monthly(i).name,...
% %         ' (',num2str(monthly(i).lat),...
% %         ', ',num2str(monthly(i).lon),')'])  

end

save('./output/monthly4.mat','monthly','soireg')

%% 6. PLOT SPATIAL VARIABILITY IN LAG CORRELATION BETWEEN ASLV AND SOI

clear all; close all;
load('./output/monthly4.mat')

h1=figure();
% Map of the Philippines
xlims=[116 128];    % longitude limits (deg E)
ylims=[4 20];       % latitude limits (deg N)
m_proj('mercator','lon',xlims,'lat',ylims); % map projection
m_gshhs_i()         % Generate coastline
m_grid
hold on

% Take only a subset - those with 0 to negative lags;
%selectrows = find(soireg(:,3) <= 0 & soireg(:,3)>-6);
selectrows = [1,4,5,6,8,9,10,14,17];

x = soireg(selectrows,1);
y = soireg(selectrows,2);
z = soireg(selectrows,3);

[xi, yi] = meshgrid(...
    linspace(min(x),max(x)),...
    linspace(min(y),max(y)));

zi = griddata(x,y,z, xi,yi);
%m_contourf(xi,yi,zi)
m_pcolor(xi,yi,zi); shading flat; colorbar
h = colorbar; ylabel(h,'lag (month)')

% Save figure to file dialog box
display(sprintf('Saving figure to file...'))
[filename,pathname]=uiputfile('*.png','Save map of stations as (*.png)');
saveas(h1,[pathname filename]);
display(sprintf('Figure saved to %s%s',pathname,filename))


%% LAG REGRESSION DETRENDED ASLV WITH SST (East)

clear all; close all

load('./output/monthly3.mat')       %SLV
load('./data_params/tao_sst.mat')    %SST

sstreg = [];
buoy = east;

for i=1:length(monthly)
    i
    t_obs = monthly(i).time;
    ind = find(t_obs >= min(buoy.time) & t_obs <= max(buoy.time));    
    [m,locb] = ismember(t_obs(ind),buoy.time);
    x = monthly(i).aslv_det2(ind);
    y = naninterp(buoy.sst(locb));
    %y = buoy.sst(locb);
    x = naninterp(monthly(i).aslv_det2(ind));
    [z,lag] = conn(x,y,10); 
    % find out at which lag correlation is highest
    [mval,mind] = max(abs(z));
    optlag = lag(mind);
    % save values to matrix
    sstreg = vertcat(sstreg,[monthly(i).lon,monthly(i).lat,...
        optlag,z(mind)]);
       
% %     figure()
% %     plot(lag,z)
% %     xlabel('lag (months)')
% %     ylabel('correlation')
% %     % note: conn(x,y,n)
% %     % x leads for + lag; y leads for - lag
% %     title([monthly(i).name,' opt lag: ',num2str(optlag),...
% %         '; %cor: ',num2str(z(mind))]);   
% %     
% %     figure()
% %     plot(t_obs,monthly(i).aslv_det2,'k-'); hold on
% %     plot(buoy.time,buoy.sst,'r-')
% %     datetick('x','mmm-yy','keeplimits','keepticks')
% %     title(['Station: ',monthly(i).name,...
% %         ' (',num2str(monthly(i).lat),...
% %         ', ',num2str(monthly(i).lon),')'])  
% % 
end
% % 
% % save('./output/monthly5.mat','monthly','sstreg')



%%  LAG REGRESSION DETRENDED ASLV WITH PDO

% % clear all; close all;
% % load('./output/monthly3.mat')   %observations
% % load('./data_params/pdo.mat')        %PDO
% % 
% % pdoreg = [];
% % for i=1:length(monthly)
% %     t_obs = monthly(i).time;
% %     ind = find(t_obs >= min(pdo.time) & t_obs <= max(pdo.time));    
% %     [m,locb] = ismember(t_obs(ind),pdo.time);
% %     x = monthly(i).aslv_det2(ind);
% %     y = pdo.value(locb);
% %     x = naninterp(monthly(i).aslv_det2(ind));
% %     [z,lag] = conn(x,y,20); 
% %     % find out at which lag correlation is highest
% %     [mval,mind] = max(abs(z));
% %     optlag = lag(mind);
% %     % save values to matrix
% %     pdoreg = vertcat(pdoreg,[monthly(i).lon,monthly(i).lat,optlag,mval]);
% %        
% %     figure()
% %     plot(lag,z)
% %     xlabel('lag (months)')
% %     ylabel('correlation')
% %     % note: conn(x,y,n)
% %     % x leads for + lag; y leads for - lag
% %     title([monthly(i).name,' opt lag: ',num2str(optlag),...
% %         '; %cor: ',num2str(mval)]);   
% %     
% % % %     figure()
% % % %     plot(t_obs,monthly(i).aslv_det2,'k-'); hold on
% % % %     plot(soi.time,soi.Value,'r-')
% % % %     datetick('x','mmm-yy','keeplimits','keepticks')
% % % %     title(['Station: ',monthly(i).name,...
% % % %         ' (',num2str(monthly(i).lat),...
% % % %         ', ',num2str(monthly(i).lon),')'])  
% % 
% % end
% % 
% % save('./output/monthly5.mat','monthly','pdoreg')

%% 7. PLOT SPATIAL VARIABILITY IN LAG CORRELATION BETWEEN ASLV AND PDO

% % % clear all; close all;
% % % load('./output/monthly5.mat')
% % % 
% % % h1=figure();
% % % % Map of the Philippines
% % % xlims=[116 128];    % longitude limits (deg E)
% % % ylims=[4 20];       % latitude limits (deg N)
% % % m_proj('mercator','lon',xlims,'lat',ylims); % map projection
% % % m_gshhs_i()         % Generate coastline
% % % m_grid
% % % hold on
% % % 
% % % % Take only a subset - those with 0 to negative lags;
% % % %selectrows = find(soireg(:,3) <= 0 & soireg(:,3)>-6);
% % % selectrows = [1,4,5,6,8,9,10,14,17];
% % % 
% % % x = pdoreg(selectrows,1);
% % % y = pdoreg(selectrows,2);
% % % z = pdoreg(selectrows,3);
% % % 
% % % [xi, yi] = meshgrid(...
% % %     linspace(min(x),max(x)),...
% % %     linspace(min(y),max(y)));
% % % 
% % % zi = griddata(x,y,z, xi,yi);
% % % %m_contourf(xi,yi,zi)
% % % m_pcolor(xi,yi,zi); shading flat; colorbar
% % % % 
% % % % % Save figure to file dialog box
% % % % display(sprintf('Saving figure to file...'))
% % % % [filename,pathname]=uiputfile('*.png','Save map of stations as (*.png)');
% % % % saveas(h1,[pathname filename]);
% % % % display(sprintf('Figure saved to %s%s',pathname,filename))