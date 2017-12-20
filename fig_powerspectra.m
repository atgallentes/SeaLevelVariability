function [ output_args ] = fig_powerspectra( avec, tvec, Y )
%UNTITLED3 Summary of this function goes here
%   Input
%       avec - signal in time domain, hourly, with gaps or irregularly
%       spaced
%       tvec - time points for observations in avec
%       Y - latitude

% % figure()
[P,F] = plomb(avec,tvec);
% % plot(F,P,'k-')
y = P;

ftide=tide_freqs(char('M2','K1','S2','O1'));  % frequencies of constituents, per hour
ftide=ftide*24;                 % convert to cycles per day
f=coriolis(mean(Y(:)));         % returns degrees
fc=f*3600/(2*pi);               % rad/s to rad/hr
fc=fc*24;

figure1 = figure;

fontsize = 10;
offy = 1000000000;  %offset for y axis max

ax1 = axes('Parent',figure1,'YScale','log','XScale','log');
box(ax1,'off');
hold(ax1,'on');

% loglog(fscale,y,'Parent',ax1,'Color',[0 0 0],'LineWidth',1)
xlabel('frequency (cpd)')
%xlim([1e-03 F(end)])
xlim([1e-03 1e01])
ylim([1 max(y)+2*offy])

linecolor = [0.5 0.5 0.5];

line([ftide(2) ftide(2)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(ftide(2),max(y)+offy,'K_1','HorizontalAlignment','Left',...
    'VerticalAlignment','Top','fontsize',fontsize)

line([ftide(4) ftide(4)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(ftide(4),max(y)+offy,'O_1','HorizontalAlignment','Right',...
    'VerticalAlignment','Top','fontsize',fontsize)

line([fc fc],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(fc,max(y)+offy,'f','HorizontalAlignment','Left',...
    'VerticalAlignment','Top','fontsize',fontsize)

line([ftide(3) ftide(3)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(ftide(3),max(y)+offy,'S_2','HorizontalAlignment','Left',...
    'VerticalAlignment','Top','fontsize',fontsize)

line([ftide(1) ftide(1)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(ftide(1),max(y)+offy,'M_2','HorizontalAlignment','Right',...
    'VerticalAlignment','Top','fontsize',fontsize)

% line([(3/2)*ftide(1) (3/2)*ftide(1)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
% text((3/2)*ftide(1),max(y)+offy,'M_3','HorizontalAlignment','Left','fontsize',fontsize)
% 
% line([2*ftide(1) 2*ftide(1)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
% text(2*ftide(1),max(y)+offy,'M_4','HorizontalAlignment','Left','fontsize',fontsize)
% 
% line([3*ftide(1) 3*ftide(1)],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
% text(3*ftide(1),max(y)+offy,'M_6','HorizontalAlignment','Left','fontsize',fontsize)

line([1/365 1/365],[min(y) max(y)+offy],'Parent',ax1,'LineStyle',':','Color','r','LineWidth',0.5)
%text(1/365,max(y),'annual','HorizontalAlignment','Left','fontsize',fontsize)

loglog(F,y,'Parent',ax1,'Color',[0 0 0],'LineWidth',1)

ax1_pos = get(ax1,'Position');
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
% 1/F(end)
% 1/F(1)
%xlim([1/F(end) 1/1e-03])
xlim([1/1e01 1/1e-03])
ylim([1 max(y)+2*offy])
set(ax2,'XDir','reverse','YScale','log','XScale','log')
xlabel('period (days)')

y2 = [min(y) max(y)+offy];
x2 = [1/F(end) 1/F(1)];
line(x2,y2,'Parent',ax2,'LineStyle','none','Color','r')

line([365 365],y2,'Parent',ax2,'LineStyle',':','Color',linecolor,'LineWidth',0.5)
text(365,max(y)+offy,'annual','HorizontalAlignment','Left',...
    'VerticalAlignment','Top','fontsize',fontsize)

hold off;

end

