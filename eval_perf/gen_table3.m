%{
gen_table3.m

This file calculates the welfare benefit of giving an agent access to
Case-Shiller Index Futures

Copyright A. Michael Sharifi
%}

function [ w, vfn_store, vfn_CSF_store, w_DIFF_store ] = ...
    gen_table3(col, age0, t_i1, city_id, ds, t_begin, t_end, plotFlag, p_mid  )

save('gen_table3_save');

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
%idx_w = all( [ ds(:,col.W) >= 0.0, ds(:,col.W) <= 20.0 ], 2 );

%%
age0_store = [ 30, 45, 60 ];

%%
if (plotFlag == 1) 
    figure; hold on;
end

%%
for t = t_begin:t_end
    
    for i_age = 1:3
        age0 = age0_store(i_age);
        
        idx1 = all( [ ds(:,col.age0) == age0, ...       % age condition
            ds(:,col.city_id) == city_id, ...           % city_id condition
            ds(:,col.year_id) == t, ...                  % year condition
            ds(:,col.hor_id) == 0, ...                   % year horizon = 0
            idx, ...
            ds(:,col.i_yi) == 1, ...                     % i_yi = 1
            ds(:,col.t_i) == t_i1, ...                    % t_i = 0
            ds(:,col.ph_i) == p_mid ], 2);                % p_i = p_mid (actual price)
        
        idx2 = all( [ ds(:,col.age0) == age0, ...       % age condition
            ds(:,col.city_id) == city_id, ...           % city_id condition
            ds(:,col.year_id) == t, ...                 % year condition
            ds(:,col.hor_id) == 0, ...                  % year horizon = 0
            idx_CSF, ...
            ds(:,col.i_yi) == 1, ...                      % i_yi = 1
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
                %w_DIFF(w_i) = w(idx_w) - w(w_i);
                % try: assume linear interpolation
                w_DIFF(w_i) = (v_CSF - vfn(w_i) ) / (vfn(idx_w) - vfn(w_i) ) ...
                    * ( w(idx_w) - w(w_i) )  ;
            end
            
            w_DIFF_store(w_i,t) = w_DIFF(w_i);
            
        end
        
        %if (i_age == 3)
        %    w_DIFF_store_tmp = w_DIFF_store(:,t);
        %    for w_i = 9:length(w_DIFF)-8
        %        w_DIFF_store(w_i,t) = mean( w_DIFF_store_tmp( w_i-8:w_i+8) );
        %    end          
        %end
        
        %guess: subplot(2, Number of plots, index in plot)
        w_idx = all( [ w>= 0.0, w <= 25.0 ], 2 );
        
        if (plotFlag == 1)
            hold on;
            t_id = t - t_begin + 1;
            %subplot(2,N_plot,t_id);
            %hold on;
            %plot(w, vfn);
            %plot(w, vfn_CSF, 'r');
            
            %subplot(2,N_plot,N_plot + t_id );
            %subplot(1, N_plot, t_id);
            subplot(2, ceil(N_plot/2), t_id);
            
            if (i_age == 1)
                plot(w(w_idx), w_DIFF(w_idx));
                title_str = sprintf('Year  %d', t_id + (t_begin + 2007 - 5 -1) );
                title( title_str );
                                
            elseif( i_age == 2)
                plot(w(w_idx), w_DIFF(w_idx), ':' );
            else
                plot(w(w_idx), w_DIFF(w_idx), '--' );
            end
            
            xlim([ 0.0, 25.0]);
           
            if ( (i_age == 3 ) && (t == t_end) )
                legend('Age 30','Age 45','Age 60','location','bestoutside');    
            end
            
        end
    end
    
end

if (plotFlag == 1)    
    %suptitle('fooo');
    
    %set(gcf,'NextPlot','add');
    %axes;
    %h = title('Equivalent Wealth Gain');
    %set(gca,'Visible','off');
    %set(h,'Visible','on');
end


end




%{
city_id, year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, w_i,
W, C, B, X, CSFp, CSFn, t_i2, V
%}