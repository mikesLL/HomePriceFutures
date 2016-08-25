%{
view_pricepath.m
Script loads simulated home price paths from 'price_results', plots them,
and saves the forecast as figures in 'figures'


inputs:
city_id: city of interest
yr_id_beg: first year of interest
yr_id_end: last year of interest
vis_id: set = 1 to display figures, =0 otherwise

Copyright A. Michael Sharifi, 2016
%}

function view_pricepath( city_id, yr_id_beg, yr_id_end, vis_id )

city_str_store = {'sd','sf','lax'};
city_str = char(city_str_store(city_id + 1));  % switch from C++ to MATLAB index
years = (2003:2013);

for yr_id= yr_id_beg:yr_id_end
    disp(yr_id);
    gen_plot(years, yr_id, city_str, vis_id);
end

end

function gen_plot(years, yr_id, city_str, vis_id)

addpath('figures');
addpath('price_results');

file_str = sprintf('%syr%dp1_file.csv', city_str, yr_id); %filename of interest

ds = xlsread(file_str);                                   %load file
T_h_max = size(ds, 2);                                    %find max horizon
h_step = 0:1:(T_h_max-1);
ds_quant = zeros(5,T_h_max);                              

for i=1:T_h_max
    ds(:,i) = max( ds(:,i), 0.0);
    ds_quant(:,i) = quantile( ds(:,i) ,[ .05 .25 .5 0.75 .95] );  % compute quantiles from sample
end

h = figure;
if vis_id == 0
    set(h,'Visible','off','Position',[0 0 800 400]);
else
    set(h,'Position',[0 0 800 400]);
end

% generate plot
hold on;
plot(h_step, ds_quant(1,:));
plot(h_step, ds_quant(2,:),':');
plot(h_step, ds_quant(3,:));
plot(h_step, ds_quant(4,:),':');
plot(h_step, ds_quant(5,:));
legend('.05 quantile', '.25 quantile', 'median', '.75 quantile', '.95 quantile',...
    'Location','southoutside','Orientation','horizontal');
ylabel('Home Price Forecast');
xlabel('steps ahead');

% store plot in figures directory
str1 = sprintf('figures/forecast_%s_%d',  char(lower(city_str)), years(yr_id) );
print(str1,'-depsc','-tiff');
print(str1,'-dpng'); 

end





