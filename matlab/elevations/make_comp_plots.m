function [] = make_comp_plots(ctdData, ncomData)
% Generates a series of plots comparing ncom and observation data;
% Salinity data is smooth using a moving average with a span of 50
%
% Input:
%    ctdData - ctdData structure created by get_ctd_data
%    ncomData - ndomData structure created by get_ncom_match_ctd
% Output:
%    An image comparing temperature and salinity differences at all ctd casts 
%    sites saved locally.
%
% lopezj - 01/12/2012
%

if ~exist('./plots')
	mkdir('plots')
end

for cast = 1:size(ctdData,2)
	figH = figure('visible','off','color','w');

	% Get index
	idx = find(isnan(ncomData(cast).data(:,1)),1)-1; 

	% Salinity vs. depth
	plot(smooth(ctdData(cast).data(:,2),50),-smooth(ctdData(cast).data(:,1),50), ...
	            'color','blue','linewidth',1.5,'marker','x'); hold on;
	plot(ncomData(cast).data(1:idx,2),ncomData(cast).data(1:idx,1),'color','red', ...
	     'linewidth', 1.5);
	set(gca,'fontsize',14);

	% Annotate graph
	t = sprintf('Lat: %f Long: %f Day: %s', ctdData(cast).lat, ctdData(cast).long, ...
	            datestr(ctdData(cast).time, 'mm/dd/yy'));
	title(t);
	xlabel('Salinity psu');
	ylabel('Depth m');
	legend('CTD', 'NCOM', 'Location', 'SouthWest');

	% Save the plot
	fn = sprintf('./plots/%s-%03d.png', 'ctd_ncom', cast);  
    iw = 1024; 
    ih = 800; 
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 iw/100.0 ih/100.0])
    print('-dpng', fn, '-r100');

    close(figH);
end
