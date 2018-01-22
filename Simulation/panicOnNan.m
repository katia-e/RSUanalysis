function panicOnNan(val)
for i=1:length(val)
    if isnan(val(i))
        error('ERR: wrong value')
    end
end

