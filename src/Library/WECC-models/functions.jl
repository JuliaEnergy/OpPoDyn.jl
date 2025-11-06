function deadband(x,low,high)
    ifelse(x < low, x-low, ifelse(x>high, x-high, 0))
end

function limiter(x, lim_low, lim_high)
    #ifelse(x<lim_low, lim_low, ifelse(x>lim_high, lim_high, x))
    clamp(x, lim_low, lim_high)
end

function lowlimit(x, lim_low)
    #ifelse(x<lim_low, lim_low,  x)
    max(x, lim_low)
end

function uplimit(x, lim_high)
    #ifelse(x>lim_high, lim_high, x)
    min(x, lim_high)
end

function LVPLogic(in, zero, Brkpt, gain)
    ifelse(in<zero, 0, ifelse(in>Brkpt, gain, gain * (in-zero)/(Brkpt-zero)))
end

function VDL(V, v1, v2, v3, v4, i1, i2, i3, i4)
    ifelse(V<v1, i1+(i2-i1)/(v2-v1)*(V-v1), ifelse(V<v2, i1+(i2-i1)/(v2-v1)*(V-v1), ifelse(V<v3, i2+(i3-i2)/(v3-v2)*(V-v2), ifelse(V<v4, i3+(i4-i3)/(v4-v3)*(V-v3), i4+(i4-i3)/(v4-v3)*(V-v4)))))
end #linear interpolation; extrapolation last two points