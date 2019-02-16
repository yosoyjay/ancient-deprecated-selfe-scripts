function []= discharge(filter)
%Script plots Bonneville Dam Discharge vs # Hours Filter Criteria Satisfied
%
%Filter options: (salmon, chinook_a, chinook_s, coho_a, coho_s, sockeye_a,
%sockeye_s, steelhead_a, steelhead_s, lamprey_S5D4, lamprey_S5D16,
%lamprey_S10D4, lamprey_S10D16, lamprey_S35D4, lamprey_S35D16).
clf
%Discharge%
load('/home/users/flierm/scripts/discharge.mat', 'june_discharge', 'september_discharge');
numDays= 14;
daily_discharge=zeros(numDays,1);
total_daily_discharge=zeros(numDays*2,1);
%Filter Hours%
directory= ('/home/users/flierm/2010');
months= {'Jun/06-13-2010_14_days' 'Jun/06-13-2010_14_days_slr' 'Sep/09-12-2010_14_days' 'Sep/09-12-2010_14_days_slr'};
startDate= {164 164 255 255};
endDate= {177 177 268 268};
month_discharge= {june_discharge june_discharge september_discharge september_discharge};
dateStructure= struct ('filePath', months, 'startDate', startDate, 'endDate', endDate, 'discharge', month_discharge);
elem_area= zeros(5,numDays*2);
col_index= 1;
B=zeros(5,1);
B_index= 1;
titles= {'Bonneville Discharge vs Hours/Day Filter Criteria Satisfied: June 2010 ' 'Bonneville Discharge vs Hours/Day Filter Criteria Satisfied: June 2010 (1.2m SLR)' 'Bonneville Discharge vs Hours/Day Filter Criteria Satisfied: September 2010' 'Bonneville Discharge vs Hours/Day Filter Criteria Satisfied: September 2010 (1.2m SLR)'};
figure_directory= sprintf('%s/discharge/%s', directory, filter);
outDIR= exist (figure_directory);
if outDIR~= 7
    mkdir (figure_directory);
end
filenames= {'June_discharge.jpg' 'June_SLR_discharge.jpg' 'September_discharge.jpg' 'September_SLR_discharge.jpg'};

for allmonths=1:2
    col_index=1;
    for day= 1:numDays
		%Discharge%
		a=24*(day-1)+1;
		b=a+23;
		integrated_discharge= cumtrapz(dateStructure(allmonths).discharge(a:b,2)); 
		daily_discharge= sum(integrated_discharge);
		total_daily_discharge(day,1)= daily_discharge;
		%Filter Hours%
		julian_day=dateStructure(allmonths).startDate+day-1;
		datafile= sprintf('%s/%s/post/daily_hab_opp/2010-%i/process/%s/regional_habitat.ascii', directory, months{allmonths}, julian_day, filter);
		%disp(datafile);
		A= importdata(datafile, ' ', 2);
		for k=[1, 4, 7, 10, 11]; %In region_habitat.ascii, 1=elev, 4=avg sal, 7=avg vel, 10=avg temp, 11=avg comb
			B(B_index, :)=A.data(1,k);
			B_index=B_index+1;
		end
		B_index= 1;
		elem_area(:,col_index)= B;
		col_index=col_index+1;
    end
    for day= 1:numDays
		%Discharge%
		a=24*(day-1)+1;
		b=a+23;
		integrated_discharge= cumtrapz(dateStructure(allmonths+2).discharge(a:b,2)); 
		daily_discharge= sum(integrated_discharge);
		total_daily_discharge(day+14,1)= daily_discharge;
		%Filter Hours%
		julian_day=dateStructure(allmonths+2).startDate+day-1;
		datafile= sprintf('%s/%s/post/daily_hab_opp/2010-%i/process/%s/regional_habitat.ascii', directory, months{allmonths+2}, julian_day, filter);
		%disp(datafile);
		A= importdata(datafile, ' ', 2);
		for k=[1, 4, 7, 10, 11]; %In region_habitat.ascii, 1=elev, 4=avg sal, 7=avg vel, 10=avg temp, 11=avg comb
			B(B_index, :)=A.data(1,k);
			B_index=B_index+1;
		end
		B_index= 1;
		elem_area(:,col_index)= B;
		col_index=col_index+1;
    end

	tdd=(total_daily_discharge)';

	%old ts.m plot #hrs of habitat opportunity
	plot(tdd, elem_area(1,:), 'r.', 'MarkerSize', 20); 
	xlabel('Bonneville Dam Discharge (cubic meters/day)', 'FontSize', 15);
	ylim([0 25]);
	ylabel('Hours/Day Filter Criteria Satisfied', 'FontSize', 15);
	title(titles{allmonths}, 'FontSize', 15);
	hold on
	plot(tdd, elem_area(2,:), 'g.', 'MarkerSize', 20);
	hold on
	plot(tdd, elem_area(3,:), 'c.', 'MarkerSize', 20);
	hold on
	plot(tdd, elem_area(4,:), 'b.', 'MarkerSize', 20);
	hold on
	plot(tdd, elem_area(5,:), 'm.', 'MarkerSize', 20);
	hold on
	hl= legend('Elevation', 'Salinity', 'Velocity', 'Temperature', 'Combined', 'Location', 'EastOutside');
	set(get(hl, 'title'), 'string', 'Filter Criteria', 'FontSize', 15);
	%save graphs
	discharge_graphs= sprintf('%s/%s', figure_directory, filenames{allmonths});
    blah=sprintf('blah_%s', allmonths);
	saveas(gcf, blah);
	clf;
end


