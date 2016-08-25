%{
Given policies xpol and returns ret, evaluate household's wealth

Copyright A. Michael Sharifi, 2016
%}

function [ wealth1 ] = eval_wealth( xpol, gret, y_inc2, csfLev )


wealth1 = ...
    gret.bonds .* xpol.bonds + gret.equities .* xpol.equities + ...
    gret.csf .* csfLev .* (xpol.csfP - xpol.csfN) + xpol.csfP + xpol.csfN + ...
    gret.ph .* ten_w(xpol.ten_i2) + y_inc2;

idx_rent = (xpol.ten_i2 <= 0);

wealth1(idx_rent) = wealth1(idx_rent) - gret.rent(idx_rent);


end


%y_inc2 = .8;
%csfLev = .1031*50.0;
%xpol = xpol_11;

%wealth1 = zeros(length(xpol),1);

%{
from code, next period wealth looks something like this:
w2 = rb*x[1] + ret[x_i2] * x[2] +  csfLev * csf_net2[ph_i2] * (x[3] - x[4] ) + 
					x[3] + x[4] +					
					(double) t_i2 * (*vf2).ph_grid[ph_i2] - (double) y_i2*.2;

need to look up:
1. csfLev: comes from calibration
2. csf_net2: equals realized price minus fut price
3. y_i2: comes from calibration
%}


%{
ret = dataset;
ret.year_id = xpol.year_id;
ret.year_act = xpol.year_act;
ret.bonds = 1.05*ones(length(ret),1);
ret.equities = 1.08*ones(length(ret),1);
%}