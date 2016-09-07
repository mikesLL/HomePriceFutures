%{
eval_perf.m
This script generates summary results and tables using data from directory
first_results. first_results includes value functions and strategies from an
investor's first horizon
High priorities include tracking how tenure varies over time, city,
and wealth and how tenure are affected by access to home price index
futures

set:
%city_id = 0 for San Diego, = 1 for San Francisco, = 2 for Los Angeles
%t_use = 4;   % =4 for 2006, =11 for year 2013
reload_data


Copyright A. Michael Sharifi 2016
%}

function eval_perf( city_id, reload_data )

param.city_str_store = {'sd','sfr','lax','bos','chi','den','mia','nym'};
param.age0_store = [30, 45, 60];
param.i_yi_mid = 1;              % index: middle income state index
param.hor0 = 0;                  % index: initial horizon 
param.p_mid = 4;                 % index: initial price 
 
city_id = 0;  % set to san diego for now
%reload_data = 0; % yes, reload data
addpath('../figures');
%addpath('../first_results');

%%
reload_data = 0;
city_str = 'sd';
%city_str = 'sf';
%city_str = 'lax';
%city_str = 'bos';
%city_str = 'chi';
%city_str = 'den';
%city_str = 'mia';
%city_str = 'nym';
ds_load = load_ds( reload_data, city_str );   % loads data from initial period value function
save('ds_load_sd','ds_load');
%save('ds_load_sf','ds_load');
%save('ds_load_lax','ds_load');
%save('ds_load_bos','ds_load');
%save('ds_load_chi','ds_load');
%save('ds_load_den','ds_load');
%save('ds_load_mia','ds_load');
%save('ds_load_nym','ds_load');

%save('eval_perf_save');

%% load column indices
col.city_id = 1;      % city index
col.age0 = 2;         % initial age
col.year_id = 3;      % beginning year
col.hor_id = 4;       % horizon
col.rho = 5; 
col.gamma = 6;
col.csfLev = 7;
col.w_n = 8;
col.t_i = 9;
col.ph_i = 10;
col.i_rent = 11;
col.i_yi = 12;
col.w_i = 13;
col.W = 14;
col.C = 15;
col.B = 16;
col.X = 17;
col.csfP = 18;
col.csfN = 19;
col.t_i2 = 20;
col.V = 21;

%% initialization
y_inc2 = .56;
csfLev = 18.8; %.1031*50.0;
city_str = char(param.city_str_store(city_id+1));
years = (2007:2013);

t_begin = 5;             % beginning time period
t_end = 11;              % ending time period
rho = 1;
w_n = 400;
p_mid = 4;               % = 2 for ph_n = 5, =4 for ph_n  = 9
N_w = 400; %2000;

idx_use = all([ (ds_load(:,col.city_id) == city_id), (ds_load(:,col.rho) == rho ), ...
    (ds_load(:,col.w_n) == w_n) ], 2);
ds = ds_load(idx_use,:);                      % ds is the relevant smaller dataset

%%
age0 = 30;            % age level of interest
W1 = .4;              % three wealth levels of interest
W2 = 1.6;
W3 = 4.0;
t_i1 = 0;             % inital condition: renter

%% table of policies over time w/o CSF  (three wealth levels)
csfFlag = 0;
xpol_01 = find_pol(col, age0, W1, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );
xpol_02 = find_pol(col, age0, W2, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );
xpol_03 = find_pol(col, age0, W3, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );

%% table of policies over time w/CSF (three wealth levels)
csfFlag = 1;
xpol_11 = find_pol(col, age0, W1, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );
xpol_12 = find_pol(col, age0, W2, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );
xpol_13 = find_pol(col, age0, W3, t_i1, csfFlag, ds, t_begin, t_end, p_mid  );

%% gret: gross return
gret = dataset;
gret.year_id = xpol_11.year_id;
gret.year_act = xpol_11.year_act;
gret.bonds = 1.05*ones(length(gret),1);
gret.equities = 1.08*ones(length(gret),1);
gret.csf = 0.0*ones(length(gret),1);  % TODO: this should come from home price fut and from home price realization
gret.rent = 0.0*ones(length(gret),1);
gret.ph = 5.0*ones(length(gret),1);

gret.equities(t_begin:t_end) = [ 0.0573, .1154, .1834, -.034, .0114, .1311, .1264, .1606  ];
gret.csf(t_begin:t_end) = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ];
gret.rent(t_begin:t_end) = [ .2, .2, .2, .2, .2, .2, .2, .2 ];

%gret.csf_sd(t_begin:t_end) = [ -0.42138, -1.48676, 0.357776, 0.151061, -0.31007, 0.349826, 0.850712, 0.0159 ];
%gret.rent_sd(t_begin:t_end) = [ .227325, .23815, .248975];

read_in = xlsread('../lax_read_in.csv');
loc_begin = find( read_in(:,1) == t_begin, 1, 'first');
loc_end = find( read_in(:,1) == t_end, 1, 'first');

gret.csf(t_begin:t_end) = read_in(loc_begin+1:loc_end+1,4) - read_in(loc_begin:loc_end,6);
%gret.rent(t_begin:t_end) = [ .227325, .23815, .248975];

%%
wealth_02 = eval_wealth( xpol_02, gret, y_inc2, csfLev );
wealth_12 = eval_wealth( xpol_12, gret, y_inc2, csfLev );

xpol_02.wealth = wealth_02;
xpol_12.wealth = wealth_12;
%%
figure;
hold on;
plot(wealth_02(t_begin:t_end));
plot(wealth_12(t_begin:t_end),'r');

%%
%{
%% plot tenure over time for W1
figure; hold on;
plot(ten_i2(1,:));
plot(ten_i2_CSF(1,:)', 'r');

%% plot tenure over time for W2
figure; hold on;
plot(ten_i2(2,:));
plot(ten_i2_CSF(2,:), 'r');

%% plot tenure over time for W3
figure; hold on;
plot(ten_i2(3,:));
plot(ten_i2_CSF(3,:), 'r');
%}
%% value funcion and wealth compensation for adding CSF
% vfn_store, vfn_CSF_store: each col is a value fn for one time period

save('eval_perf2');

%% gen_table3: plots estimated wealth gains by year and for all age groups
plotFlag = 1;
%[w, vfn_store, vfn_CSF_store, w_DIFF_store ] = ...
%    gen_table3(col, age0, t_i1, city_id, ds, t_begin, t_end, plotFlag, p_mid);
[w, vfn_store, vfn_CSF_store, w_DIFF_store ] = ...
    gen_fig3(param, col, t_i1, city_id, ds, t_begin, t_end, plotFlag, p_mid);

%% lowest level of wealth at which a renter decides to become an owner
[  ten_i2, ten_i2_CSF, csf, B, B_CSF, X, X_CSF ] = ...
    gen_table4(col, age0, t_i1, city_id, ds, t_begin, t_end, p_mid  );


% note: city_id = 2 so currently looking at Los Angeles
%t_use = 4;   % =4 for 2006, =11 for year 2013
%t_use = t_use_in;


%% plot tenure choices by wealth in 2013
% tenure.png
t_use = 8;
figure('Position', [0 0 1200 400]);
h1(1) = subplot(1,2,1); 
figure; hold on;
plot(w, ten_i2(:,t_use));       %plot tenure w/o CSF by wealth
plot(w, ten_i2_CSF(:,t_use),'--');
title('Tenure w/o CSF');
xlabel('Wealth');
ylim([0, 3]);

h1(2) = subplot(1,2,2);          %plot tenure w/ CSF by wealth
plot(w, ten_i2_CSF(:,t_use));
title('Tenure w/ CSF');
xlabel('Wealth');
ylim([ 0, 3]);
%linkaxes(h1,'y', 'x');

%% Find wealth levels at which tenure switches
[idx1, ~] = find( ten_i2(:,t_use) == 1, 1, 'first');  % switch from renter to owner
[idx2, ~] = find( ten_i2(:,t_use) == 2, 1, 'first');  % switch from starter to full-size
[idx1_CSF, ~] = find( ten_i2_CSF(:,t_use) == 1, 1, 'first');  % switch from renter to owner
[idx2_CSF, ~] = find( ten_i2_CSF(:,t_use) == 2, 1, 'first');  % switch from starter to full-size

%% plot asset allocation holdings by wealth in 2013
% asset_allocation.png
figure('Position', [0 0 1200 400]);
h2(1) = subplot(1,2,1);
hold on;
plot(w(1:N_w), B(1:N_w,t_use));
plot(w(1:N_w), X(1:N_w,t_use),'.');
read_in = [ B(:,t_use) ; X(:,t_use); B_CSF(:,t_use); X_CSF(:,t_use) ];
y_upper = max( read_in(:) ) + .1*abs(max( read_in(:) ));
y_lower = min( read_in(:) ) - .1*abs(max( read_in(:) ));
if ~isempty(w(idx1))
    line([w(idx1) w(idx1)],[y_lower, y_upper],'LineStyle','--');
end

if ~isempty(w(idx2))
    line([w(idx2) w(idx2)],[y_lower, y_upper],'LineStyle','--');
end
title('Asset Allocation');
xlabel('Wealth');

h2(2) = subplot(1,2,2);
hold on;
plot(w(1:N_w), B_CSF(1:N_w,t_use));
plot(w(1:N_w), X_CSF(1:N_w,t_use),'.');
if ~isempty(w(idx1_CSF))
    line([w(idx1_CSF) w(idx1_CSF)],[y_lower, y_upper],'LineStyle','--');
end

if ~isempty(w(idx2_CSF))
    line([w(idx2_CSF) w(idx2_CSF)],[y_lower, y_upper],'LineStyle','--');
end
title('Asset Allocation');
xlabel('Wealth');
linkaxes(h2,'y');


%% CSF plots here
%print(strcat('figs/negeq_',char(lower(city_str))),'-depsc','-tiff'); % Neg Eq homeowners
%print(strcat('figs/negeq_',char(lower(city_str))),'-dpng'); 
 

%% plot CSF positions by wealth in 2013
% csf_position.png
figure('Position', [0 0 600 400]);
plot(w(1:N_w), csf(1:N_w,t_use));
yL = get(gca,'YLim');
if ~isempty(w(idx1_CSF))
    line([w(idx1_CSF) w(idx1_CSF)],yL,'LineStyle','--');
end

if ~isempty(w(idx2_CSF))
    line([w(idx2_CSF) w(idx2_CSF)],yL,'LineStyle','--');
end
title('CSF Position (allocated margin)');
xlabel('Wealth');

str1 = sprintf('figures/csf_position_%s_%d',city_str,years(t_use));
print(str1,'-depsc','-tiff'); % csf position
print(str1,'-dpng'); % csf position

%% plot Value Function with and without CSF
% value_fn.png
figure('Position', [0 0 1200 400]);
subplot(1,2,1);
hold on;
plot(w, vfn_store(:,t_use));
plot(w, vfn_CSF_store(:,t_use), '--');
title('Value Functions');
legend('w/o CSF', 'w/ CSF', 'Location', 'southeast' );
xlabel('Wealth');

subplot(1,2,2);
plot(w, w_DIFF_store(:,t_use) );
yL = get(gca,'YLim');
if ~isempty(w(idx1_CSF))
    line([w(idx1_CSF) w(idx1_CSF)],yL,'LineStyle','--');
end

if ~isempty(w(idx2_CSF))
    line([w(idx2_CSF) w(idx2_CSF)],yL,'LineStyle','--');
end
title('Equivalent Wealth Gain');
xlabel('Wealth');


%{
city_id, year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, w_i,
W, C, B, X, CSFp, CSFn, t_i2, V
%}

%% Table: wealth levels at which tenure switches
w_switch = zeros(2,t_end);
for t = t_begin:t_end
    [idx1, ~] = find( ten_i2(:,t) == 1, 1, 'first');  % switch from renter to owner
    [idx1_CSF, ~] = find( ten_i2_CSF(:,t) == 1, 1, 'first');  % switch from renter to owner
    
    if ~isempty(w(idx1))
        w_switch(1, t) = w(idx1);
    end
    
    if ~isempty(w(idx1_CSF))
        w_switch(2, t) = w(idx1_CSF);
    end
end

yr = (2006:2013);
fore = [-.1536, .0078, -.1474, .0677, .1244, .079, .1241, .1093];

figure;
plot(yr, w_switch(:,4:11)');

figure;
plot(yr, w_switch(2,4:11) ./ w_switch(1,4:11));

figure;
plot(yr, fore);


end





