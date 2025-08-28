function deadband(x,low,high)
    ifelse(x < low, x-low, ifelse(x>high, x-high, 0))
end

function limiter(x, lim_low, lim_high)
    ifelse(x<lim_low, lim_low, ifelse(x>lim_high, lim_high, x))
end

function lowlimit(x, lim_low)
    ifelse(x<lim_low, lim_low,  x)
end

function uplimit(x, lim_high)
    ifelse(x>lim_high, lim_high, x)
end

function LVPLogic(in, zero, Brkpt, gain)
    ifelse(in<=Brkpt, ifelse(in>=zero, gain * (in-zero)/(Brkpt-zero), 0), 0)
end
