%{
gen_table1.m

This file plots the value function with wealth on the x axis given current 
state variables and access to home price futures

Copyright A. Michael Sharifi
%}

function [ w, vfn ] = gen_table1(col, t_i1, csfFlag, ds, t_begin, t_end, p_mid  )

if csfFlag == 1
   idx_CSF = ( ds(:,col.csfLev) > 0.0 );
else
    idx_CSF = (ds(:,col.csfLev) <= 0.0 );
end

w_n = ds(1,col.w_n);
w = zeros(w_n, t_end);
vfn = zeros(w_n, t_end);

%%
for t = t_begin:t_end
    idx2 = all( [ ds(:,col.year_id) == t, ...        % year condition
                  ds(:,col.hor_id) == 0, ...        % year horizon = 0
                  idx_CSF, ...
                  ds(:,col.t_i) == t_i1, ...     % t_i = 0
                  ds(:,col.ph_i) == p_mid ], 2);      % p_i = 2 (actual price)
    w(:,t) = ds(idx2, col.W)';
    vfn(:,t) = ds(idx2,col.V)';
end

end
