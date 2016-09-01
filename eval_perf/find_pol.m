%{
find_pol.m
Script finds optimal tenure / home size as a function of current 
state variables and access to home price futures

Copyright A. Michael Sharifi, 2016

inputs:
W1: wealth level of interest
t_i1: current tenure
csfFlag: allow / disallow home price futures
ds: input dataset
%}

function [ xpol ] = find_pol(col, age0, W1, t_i1, csfFlag, ds, t_begin, t_end, p_mid  )

%save('find_pol_save');

% find wealth index W_i1 associated with W1 
ds_W = sortrows( ds(:,[col.w_i,col.W]) );
idx0 = find(ds_W(:,2) >= W1, 1, 'first');
W_i1 = ds_W(idx0,1);

xpol = dataset;
xpol.year_id = zeros(t_end,1);
xpol.year_act = zeros(t_end,1);
xpol.cons = zeros(t_end,1);        % consumption
xpol.bonds = zeros(t_end,1);       % bonds
xpol.equities = zeros(t_end,1);    % equity holdings
xpol.csfP = zeros(t_end,1);        % home price futures: long 
xpol.csfN = zeros(t_end,1);        % home price futures: short
xpol.ten_i2 = zeros(t_end,1);      % next-period tenure

if csfFlag == 1
   idx_CSF = ( ds(:,col.csfLev) > 0.0 );
else
    idx_CSF = (ds(:,col.csfLev) <= 0.0 );
end

%%
for t = t_begin:t_end
    idx2 = all( [ ds(:,col.age0) == age0, ...   % age condition
        ds(:,col.year_id) == t, ...             % year condition
        ds(:,col.hor_id) == 0, ...              % year horizon = 0
        idx_CSF, ...
        ds(:,col.t_i) == t_i1, ...              % t_i = 0
        ds(:,col.i_yi) == 1, ...              % i_yi = 1
        ds(:,col.ph_i) == p_mid, ...            % p_mid: represents actual price
        ds(:,col.w_i) == W_i1] , ...            % find w_i of interest
        2);
    
    xpol.year_id(t) = t;
    xpol.year_act(t) = t + 2002;
    xpol.cons(t) = ds(idx2,col.C);           % consumption
    xpol.bonds(t) = ds(idx2,col.B);          % bonds
    xpol.equities(t) = ds(idx2,col.X);       % equity holdings
    xpol.csfP(t) = ds(idx2,col.csfP);        % home price futures: long
    xpol.csfN(t) = ds(idx2,col.csfN);        % home price futures: short
    xpol.ten_i2(t) = ds(idx2,col.t_i2);      % next-period tenure
end

end
