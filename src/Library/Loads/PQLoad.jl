@mtkmodel PQLoad begin
    @parameters begin
        P_set, [description="Active Power demand"]
        Q_set, [description="Reactive Power demand"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        P(t), [description="Active Power"]
        Q(t), [description="Reactive Power"]
    end
    @equations begin
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
        # makes it easier for the solver then P ~ Pset and Q ~ Qset
        terminal.i_r ~  simplify(real((P_set + im*Q_set)/(terminal.u_r + im*terminal.u_i)))
        terminal.i_i ~ -simplify(imag((P_set + im*Q_set)/(terminal.u_r + im*terminal.u_i)))
    end
end

@mtkmodel VoltageDependentLoad begin
    @parameters begin
        P_set, [description="Active Power demand"]
        Q_set, [description="Reactive Power demand"]
        α_P, [description="Active Power exponent"]
        α_Q, [description="Reactive Power exponent"]
        V_n, [description="Nominal voltage (where real power equals set power)"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        Q(t), [description="Reactive Power [pu]"]
    end
    begin
        v = sqrt(terminal.u_r^2 + terminal.u_i^2)
        Pload = P_set * (v/V_n)^α_P
        Qload = Q_set * (v/V_n)^α_Q
        Sload = Pload + im*Qload
        vcomplex = terminal.u_r + im*terminal.u_i
        iout = conj(Sload/vcomplex)
    end
    @equations begin
        terminal.i_r ~ simplify(real(iout))
        terminal.i_i ~ simplify(imag(iout))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end
end

# FIXME: handle Q_Set = 0 and Q_set != 0 in the same model?
@mtkmodel ConstantYLoad begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        P_set, [description="Active Power demand [pu]"]
        Q_set, [description="Reactive Power demand [pu]"]
        V_set, [guess=1,description="Nominal voltage [pu]"]
    end
    @variables begin
        P(t), [description="Active Power [pu]"]
        Q(t), [description="Reactive Power [pu]"]
    end
    begin
        S = P_set + im*Q_set #TODO: currently not working with Q=0
        #S = Pset
        Y = conj(S)/V_set^2
        iload = Y * (terminal.u_r + im*terminal.u_i)
    end
    @equations begin
        terminal.i_r ~ simplify(real(iload))
        terminal.i_i ~ simplify(imag(iload))

        # observables
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i
    end
end

# TODO: S_b for loads?
@mtkmodel ZIPLoad begin
    @parameters begin
        P_set, [description="Active Power at operation point [pu]"]
        Q_set, [description="Reactive Power at operation point [pu]"]
        V_set, [guess=1,description="Voltage at operation point [pu]"]
        K_pZ, [description="Active power constant impedance fraction"]
        K_qZ, [description="Reactive power constant impedance fraction"]
        K_pI, [description="Active power constant current fraction"]
        K_qI, [description="Reactive power constant current fraction"]
        K_pC=1-K_pZ-K_pI, [description="Active power constant power fraction"]
        K_qC=1-K_qZ-K_qI, [description="Reactive power constant power fraction"]
    end
    @components begin
        terminal = Terminal()
    end
    @variables begin
        V_rel(t), [description="Relative voltage magnitude"]
        P(t), [description="Active Power"]
        Q(t), [description="Reactive Power"]
    end
    @equations begin
        V_rel ~ sqrt(terminal.u_r^2 + terminal.u_i^2)/V_set
        P ~ P_set*(K_pZ*V_rel^2 + K_pI*V_rel + K_pC)
        Q ~ Q_set*(K_qZ*V_rel^2 + K_qI*V_rel + K_qC)
        # formulate equations for i_r and i_i instead
        terminal.i_r ~  simplify(real((P + im*Q)/(terminal.u_r + im*terminal.u_i)))
        terminal.i_i ~ -simplify(imag((P + im*Q)/(terminal.u_r + im*terminal.u_i)))
    end
end
