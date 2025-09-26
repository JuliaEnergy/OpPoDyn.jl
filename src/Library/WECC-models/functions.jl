function deadband(x,low,high)
    ifelse(x < low, x-low, ifelse(x>high, x-high, 0))
end #TODO überprüfen: ist das korrekt?

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
    ifelse(in<=Brkpt, ifelse(in>=zero, gain * (in-zero)/(Brkpt-zero), 0), gain)
end

function VDL(V, v1, v2, v3, v4, i1, i2, i3, i4)
    ifelse(V<v1, i1, ifelse(V<v2, i1+(i2-i1)/(v2-v1)*(V-v1), ifelse(V<v3, i2+(i3-i2)/(v3-v2)*(V-v2), ifelse(V<v4, i3+(i4-i3)/(v4-v3)*(V-v3), i4))))
end