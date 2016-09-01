%{
gen_table2.m

This file generates optimal tenure / home size as a function of current 
state variables and access to home price futures

Copyright A. Michael Sharifi
%}

function [ ten_i2 ] = gen_table2(col, W1, t_i1, csfFlag, ds, t_begin, t_end, p_mid  )

ds_W = sortrows( ds(:,[col.w_i,col.W]) );
idx0 = find(ds_W(:,2) >= W1, 1, 'first');
W_i1 = ds_W(idx0,1);  % W_i1: w_i (wealth index) associated with W1

ten_i2 = zeros(1,t_end);

if csfFlag == 1
   idx_CSF = ( ds(:,col.csfLev) > 0.0 );
else
    idx_CSF = (ds(:,col.csfLev) <= 0.0 );
end

%%
for t = t_begin:t_end
    idx2 = all( [ ds(:,col.year_id) == t, ...         % year condition
                  ds(:,col.hor_id) == 0, ...          % year horizon = 0
                  idx_CSF, ...
                  ds(:,col.t_i) == t_i1, ...          % t_i = 0
                  ds(:,col.ph_i) == p_mid, ...        % p_i = 2 (actual price)
                  ds(:,col.w_i) == W_i1] , ...        % find w_i of interest
                  2);   
    ten_i2(t) = ds(idx2,col.t_i2);
end

end
