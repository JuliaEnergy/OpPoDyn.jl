@mtkmodel P_I_Lim begin
    @structural_parameters begin
        K_p # Gain proportional
        K_i # Gain integral
        T=1 # Time constant
        outMin #lower limit
        outMax #upper limit
        guess=0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=guess, description="Output signal", output=true]
        x(t), [description="Integral variable"]
    end
    @equations begin
        T*Dt(x) ~ K_i * in
        out ~ K_p*in + clamp(x, outMin, outMax)
    end
end

@mtkmodel P_I_Lim_freeze begin
    @structural_parameters begin
        K_p # Gain proportional
        K_i # Gain integral
        T=1 # Time constant
        outMin #lower limit
        outMax #upper limit
        freeze #1 if integrator state is frozen, 0 else
        guess=0
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t), [guess=guess, description="Output signal", output=true]
        x(t), [description="Integral variable"]
    end
    @equations begin
        T*Dt(x) ~ ifelse(freeze==1, 0, K_i * in)
        out ~ K_p*in + clamp(x, outMin, outMax)
    end
end

@mtkmodel SimpleLag_freeze begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        freeze #1 if integrator state is frozen, 0 else
        guess=0
        default=nothing
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t)=default, [guess=guess, description="Output signal", output=true]
        x(t), [description="Internal variable"]
    end
    @equations begin
        if T < 0
            error("Time constant in SimpleLag_freeze is less than zero")
        else
            Dt(x) ~ ifelse(T==0, 0, ifelse(freeze==1, 0 , (K*in - x)/T))
            out ~ ifelse(T==0, in, x)
        end
    end
end

@mtkmodel SimpleLag_2MaxLims begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        doutMax #upper limit for deviation
        outMax #upper limit
        guess=0
        default=nothing
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t)=default, [guess=guess, description="Output signal", output=true]
        f(t), [description="Internal variable"]
        x(t), [description="Internal limited variable"]
    end
    @equations begin
        if T < 0
            error("Time constant in SimpleLag_2MaxLims is less than zero")
        else
            f ~ ifelse(T==0, 0, (K*in - x)/T)
            Dt(x) ~ min(f, doutMax)
            out ~ ifelse(T==0, min(in, outMax), min(x, outMax))
        end
    end
end

@mtkmodel SimpleLag_2Lims_freeze begin
    @structural_parameters begin
        K # Gain
        T # Time constant
        freeze #1 if integrator state is frozen, 0 else
        doutMin #lower limit for deviation
        doutMax #upper limit for deviation
        outMin #lower limit
        outMax #upper limit
        guess=0
        default=nothing
    end
    @variables begin
        in(t), [description="Input signal", input=true]
        out(t)=default, [guess=guess, description="Output signal", output=true]
        f(t), [description="Internal variable"]
        x(t), [description="Internal limited variable"]
    end
    @equations begin
        if T < 0
            error("Time constant in SimpleLag_2Lims_freeze is less than zero")
        else
            f ~ ifelse(T==0, 0, ifelse(freeze==1, 0 , (K*in - x)/T))
            Dt(x) ~ clamp(f, doutMin, doutMax)
            out ~ ifelse(T==0, clamp(in, outMin, outMax), clamp(x, outMin, outMax))
        end
    end
end