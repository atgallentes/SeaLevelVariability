function [ site ] = calc_tides( aslv, ttime, lat )
%CALC_TIDES Analysis of tidal constituents using t_tide functions
% INPUT
%   aslv - sea level anomaly values, regularly sampled, with or without NaN
%   ttime - datenum vector of observation times
% OUTPUT
%   site - all tidal info; see site.readme below


tpres = aslv';  %note: tpres must be a row vector in this case

% initialize matrices:
%constituents={'M2','S2','K1','O1'};
constituents={'M2','S2','N2','K2','K1','O1','Q1','P1','M3','M4','M6'};
[m,n]=size(tpres);
p=length(constituents);
Tamp=NaN*ones(m,p);
Tampci=Tamp;
Tpha=Tamp;
Tphaci=Tamp;
Ttides=NaN*ones(m,n);
Tdata=zeros(m,1);

for i=1:m
    xseries=tpres(i,:);
    I=find(~isnan(xseries));
    Tdata(i)=length(I);
    if ~isempty(I)
        try
        [name,freq,tidecon,tides]=t_tide(xseries,'interval',1,'start time',ttime(1),'latitude',lat,'rayleigh',constituents,'synthesis',0); 
        Tamp(i,:)=tidecon(:,1)';
        Tampci(i,:)=tidecon(:,2)';
        Tpha(i,:)=tidecon(:,3)';
        Tphaci(i,:)=tidecon(:,4)';
        Ttides(i,:)=tides;
        catch
        end
    end
end

% site.constituents = {'O1','K1','M2','S2'};
% site.relcon = (Tamp(3)+Tamp(4))/(Tamp(1)+Tamp(2))
site.constituents = {'Q1','O1','P1','K1','N2','M2','S2','K2',...
    'M3','M4','M6'};
site.relcon = (Tamp(6)+Tamp(7))/(Tamp(2)+Tamp(4))
site.amp = Tamp;
site.ampci = Tampci;
site.pha = Tpha;
site.phaci = Tphaci;
site.tide = Ttides;
site.T = ttime;
site.data = Tdata;
site.readme=char('T_tide analysis of tide gauge hourly data',......
    'amp: amplitude (m/s)',......
    'ampci: 95% CI on amplitude (m/s)',......
    'pha: zonal Greenwich phase (degrees)',......
    'phaci: 95% CI on zonal Greenwich phase (degrees)',......
    'tides: tidal currents predicted (dbar)',......  
    'T: time for tidal currents (datenum format)',...
    'data: number of observations used for the harmonic analysis',......
    'Remark: nodal corrections have been applied.');
    % 
    % 'aslv: sealevel anomaly',... 
end

