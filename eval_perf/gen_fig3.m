%{
gen_fig3.m
This file calculates the welfare benefit of giving an agent access to
Case-Shiller Index Futures
Copyright A. Michael Sharifi

%}

function [ w, vfn_store, vfn_CSF_store, w_DIFF_store ] = ...
    gen_fig3(param, col, t_i1, city_id, ds, t_begin, t_end, plotFlag, p_mid  )

save('gen_fig3_save');

%%
N_plot = t_end - t_begin + 1;  % number of plots
w_n = ds(1,col.w_n);
w = zeros(w_n, 1);
vfn = zeros(w_n, 1);           % value function without CSF
vfn_CSF = zeros(w_n, 1);       % value function with CSF
w_DIFF = zeros(w_n, 1);        % difference in wealth

vfn_store = zeros(w_n, t_end);
vfn_CSF_store = zeros(w_n, t_end);
w_DIFF_store = zeros(w_n, t_end);

%%
if (plotFlag == 1) 
    h = figure; hold on;
end

%index stores current city, initial horizon, default income, default home
%price
idx0 = all( [ ds(:,col.city_id) == city_id, ...       % city_id condition
            ds(:,col.hor_id) == param.hor0, ...           % year horizon = 0
            ds(:,col.i_yi) == param.i_yi_mid, ...         % i_yi = 1
            ds(:,col.t_i) == t_i1, ...                    % t_i = 0
            ds(:,col.ph_i) == p_mid ], 2);                % p_i = p_mid (actual price)

%% cycle through years
for t = t_begin:t_end
    
    for i_age = 1:length(param.age0_store)
        age0 = param.age0_store(i_age);
        
        idx1 = all( [ idx0, ...
            ds(:,col.age0) == age0, ...                   % age condition
            ds(:,col.year_id) == t, ...                   % year condition
            ds(:,col.csfLev) <= 0.0  ], 2);               % csf Leverage
        
        idx2 = all( [ idx0, ...
            ds(:,col.age0) == age0, ...                   % age condition
            ds(:,col.year_id) == t, ...                   % year condition
            ds(:,col.csfLev) > 0.0  ], 2);               % csf Leverage
       
        w(:,1) = ds(idx1, col.W);
        vfn(:,1) = ds(idx1,col.V);
        vfn_CSF(:,1) = ds(idx2,col.V);
        vfn_store(:,t) = vfn(:,1);
        vfn_CSF_store(:,t) = vfn_CSF(:,1);
        
        % for each wealth w_i, calculate how much addl wealth must be given
        % to agent who does not have access to csf in order to have same
        % value 
        for w_i = 1:length(w_DIFF)
            v_CSF = vfn_CSF(w_i);
            [idx_w, ~ ] = find(vfn >= v_CSF, 1, 'first');
            if ( isempty(idx_w) || (idx_w - w_i <= 0 ) )
                w_DIFF(w_i) = 0.0;
            else
                % assume linear interpolation
                %w_DIFF(w_i) = (v_CSF - vfn(w_i) ) / (vfn(idx_w) - vfn(w_i) ) ...
                %    * ( w(idx_w) - w(w_i) )  ;
                
                w_DIFF(w_i) = ( w(idx_w-1) - w(w_i) ) + ...
                    (v_CSF - vfn(idx_w-1) ) / (vfn(idx_w) - vfn(idx_w-1) ) ...
                    * ( w(idx_w) - w(idx_w-1) )  ;
                
                
            end
            w_DIFF_store(w_i,t) = w_DIFF(w_i);
            
        end
        
        % restrict plot to subset of wealth 
        w_idx = all( [ w>= 0.0, w <= 10.0 ], 2 );
        
        if (plotFlag == 1)
            hold on;
            t_id = t - t_begin + 1;
            h(t_id) = subplot(2, ceil(N_plot/2), t_id);
            
            if (i_age == 1)
                plot(w(w_idx), w_DIFF(w_idx));
                title_str = sprintf('Year  %d', t_id + (t_begin + 2007 - 5 -1) );
                title( title_str );
                xlabel( 'Wealth ($100k)' );
                ylabel( 'Welfare Gain ($100k)' );
                                
            elseif( i_age == 2)
                plot(w(w_idx), w_DIFF(w_idx), ':' );
            else
                plot(w(w_idx), w_DIFF(w_idx), '--' );
            end
            
            %xlim([ 0.0, 10.0]);
           
            if ( (i_age == 3 ) && (t == t_end) )
                legend('Age 30','Age 45','Age 60','location','bestoutside');    
            end
            
        end
    end
    
end

linkaxes(h);
 ylim([0 0.6])


end
