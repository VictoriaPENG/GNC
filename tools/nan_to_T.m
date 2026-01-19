function y = nan_to_T(x, T)
%NAN_TO_T Replace NaN with T (scalar helper for post-processing metrics).
if isnan(x)
    y = T;
else
    y = x;
end
end