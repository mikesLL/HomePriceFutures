%{
gen_table3.m

This file calculates the welfare benefit of giving an agent access to
Case-Shiller Index Futures

Copyright A. Michael Sharifi
%}

function [ w, vfn_store, vfn_CSF_store, w_DIFF_store ] = ...
    gen_table3(col, t_i1, city_id, ds, t_begin, t_end, plotFlag, p_mid  )

N_plot = t_end - t_begin + 1;
w_n = ds(1,col.w_n);
w = zeros(w_n, 1);
vfn = zeros(w_n, 1);
vfn_CSF = zeros(w_n, 1);
w_DIFF = zeros(w_n, 1);

vfn_store = zeros(w_n, t_end);
vfn_CSF_store = zeros(w_n, t_end);
w_DIFF_store = zeros(w_n, t_end);

%%
idx = (ds(:,col.csfLev) <= 0.0 );
idx_CSF = (ds(:,col.csfLev) > 0.0 );     

%%
for t = t_begin:t_end
    idx1 = all( [ ds(:,col.city_id) == city_id, ...  % city_id condition
        ds(:,col.year_id) == t, ...                  % year condition
        ds(:,col.hor_id) == 0, ...                   % year horizon = 0
        idx, ...
        ds(:,col.t_i) == t_i1, ...                     % t_i = 0
        ds(:,col.ph_i) == p_mid ], 2);                 % p_i = p_mid (actual price)
    
    idx2 = all( [ ds(:,col.city_id) == city_id, ...    % city_id condition
        ds(:,col.year_id) == t, ...                    % year condition
        ds(:,col.hor_id) == 0, ...                     % year horizon = 0
        idx_CSF, ...
        ds(:,col.t_i) == t_i1, ...                    % t_i = 0
        ds(:,col.ph_i) == p_mid ], 2);                % p_i = p_mid (actual price)
    
    w(:,1) = ds(idx1, col.W);
    vfn(:,1) = ds(idx1,col.V);
    vfn_CSF(:,1) = ds(idx2,col.V);
    vfn_store(:,t) = vfn(:,1);
    vfn_CSF_store(:,t) = vfn_CSF(:,1);
    
    for w_i = 1:length(w_DIFF)
        v_CSF = vfn_CSF(w_i);
        [idx_w, ~ ] = find(vfn >= v_CSF, 1, 'first');
        if isempty(idx_w)
            w_DIFF(w_i) = 0.0;
        else
            w_DIFF(w_i) = w(idx_w) - w(w_i);
        end
        
        w_DIFF_store(w_i,t) = w_DIFF(w_i);
        
    end
   
    if (plotFlag == 1) 
        t_id = t - t_begin + 1;
        subplot(2,N_plot,t_id);
        hold on;
        plot(w, vfn);
        plot(w, vfn_CSF, 'r');
        
        subplot(2,N_plot,N_plot + t_id );
        plot(w, w_DIFF);
    end

end

end




%{
city_id, year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, w_i,
W, C, B, X, CSFp, CSFn, t_i2, V
%}