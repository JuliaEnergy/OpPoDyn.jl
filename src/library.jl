import ModelingToolkitStandardLibrary.Blocks: RealInput
import ModelingToolkit: t_nounits as t

export DQBus, DQSwing, DQLine, DQPiLine

# @connector Terminal begin
#     u_r(t), [description="d-voltage"]
#     u_i(t), [description="q-voltage"]
#     i_r(t), [description="d-current", connect=Flow]
#     i_i(t), [description="q-current", connect=Flow]
# end
@mtkmodel Terminal begin
    @variables begin
        u_r(t), [description="d-voltage"]
        u_i(t), [description="q-voltage"]
        i_r(t), [description="d-current"]
        i_i(t), [description="q-current"]
    end
end

@mtkmodel SauerPaiMachine begin
    @components begin
        terminal=Terminal()
    end
    @parameters begin
        ω_s=2π*50, [description="synchronous speed"]
        R_s, [description="stator resistance"]
        X_d, [description="d-axis synchronous reactance"]
        X′_d, [description="d-axis transient reactance"]
        X_q, [description="q-axis synchronous reactance"]
        X′_q, [description="q-axis transient reactance"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        T′_d0, [description="d-axis transient time constant"]
        T″_d0, [description="d-axis subtransient time constant"]
        T′_q0, [description="q-axis transient time constant"]
        T″_q0, [description="q-axis subtransient time constant"]
        H, [description="inertia constant"]
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        # ψ_0(t), [description="zero sequence flux linkage"]
        ψ_1d(t), [description="flux linkage assosciated with X″_d"]
        ψ_2q(t), [description="flux linkage assosciated with X″_q"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        # I_0(t), [description="zero sequence current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        # V_0(t), [description="zero sequence voltage"]
        E_fd(t), [description="field voltage"]
        E′_d(t), [description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [description="transient voltage behind transient reactance in q-axis"]
        δ(t), [description="rotor angle"]
        ω(t), [description="rotor speed"]
        T_m(t), [description="mechanical torque"]
        T_fw(t), [description="windage and friction torque"]
    end
    @equations begin
        1/ω_s * Dt(ψ_d) ~ R_s*I_d + ω/ω_s * ψ_q + V_d
        1/ω_s * Dt(ψ_q) ~ R_s*I_q - ω/ω_s * ψ_d + V_q
        # 1/ω_s * Dt(ψ_0) ~ R_s*I_0 + V_0
        T′_d0 * Dt(E′_q) ~ -E′_q - (X_d - X′_d)* (I_d - (X′_d-X″_d)/(X′_d-X_ls)^2 * (ψ_1d + (X′_d-X_ls)*I_d - E′_q)) + E_fd
        T″_d0 * Dt(ψ_1d) ~ -ψ_1d + E′_q - (X′_d - X_ls)*I_d
        T′_q0 * Dt(E′_d) ~ -E′_d - (X_q - X′_q)* (I_q - (X′_q-X″_q)/(X′_q-X_ls)^2 * (ψ_2q + (X′_q-X_ls)*I_q - E′_d))
        T″_q0 * Dt(ψ_2q) ~ -ψ_2q - E′_d - (X′_q - X_ls)*I_q
        Dt(δ) ~ ω - ω_s
        2*H/ω_s * Dt(ω) ~ T_m  - (ψ_d*I_q - ψ_q*I_d) - T_fw
        ψ_d ~ -X″_d*I_d + (X″_d - X_ls)/(X′_d - X_ls) * E′_q + (X′_d - X″_d)/(X′_d - X_ls) * ψ_1d
        ψ_q ~ -X″_q*I_q + (X″_q - X_ls)/(X′_q - X_ls) * E′_d + (X′_q - X″_q)/(X′_q - X_ls) * ψ_2q
        # ψ_0 ~ -X_ls*I_0+ V_0
    end
end

export SystemBase
@mtkmodel SystemBase begin
    @parameters begin
        SnRef = 100, [description="System base"]
        fNom = 50, [description="AC system frequency"]
        ωNom = 2 * π * fNom, [description="System angular frequency"]
        ωRef0Pu = 1, [description="Reference for system angular frequency (pu base ωNom)"]
        ω0Pu = 1, [description="System angular frequency (pu base ωNom)"]
    end
end

export DynawoMachine
@mtkmodel DynawoMachineParameterBase begin
    @parameters begin
        # General parameters of the synchronous machine
        UNom=24, [description="Nominal voltage in kV"]
        SNom=2220, [description="Nominal apparent power in MVA"]
        PNomTurb=2220, [description="Nominal active (turbine) power in MW"]
        PNomAlt=2200, [description="Nominal active (alternator) power in MW"]
        QNomAlt=sqrt(SNom * SNom - PNomAlt * PNomAlt), [description="Nominal reactive (alternator) power in Mvar"]
        ExcitationPu, [description="Choice of excitation base voltage"]
        H=3.5, [description="Kinetic constant = kinetic energy / rated power"]
        DPu=0, [description="Damping coefficient of the swing equation in pu"]

        # Transformer input parameters
        SnTfo=2220, [description="Nominal apparent power of the generator transformer in MVA"]
        UNomHV=24, [description="Nominal voltage on the network side of the transformer in kV"]
        UNomLV=24, [description="Nominal voltage on the generator side of the transformer in kV"]
        UBaseHV=24, [description="Base voltage on the network side of the transformer in kV"]
        UBaseLV=24, [description="Base voltage on the generator side of the transformer in kV"]
        RTfPu=0, [description="Resistance of the generator transformer in pu (base UBaseHV, SnTfo)"]
        XTfPu=0, [description="Reactance of the generator transformer in pu (base UBaseHV, SnTfo)"]

        # Mutual inductances saturation parameters, Shackshaft modelisation
        md=0.031, [description="Parameter for direct axis mutual inductance saturation modelling"]
        mq=0.031, [description="Parameter for quadrature axis mutual inductance saturation modelling"]
        nd=6.93, [description="Parameter for direct axis mutual inductance saturation modelling"]
        nq=6.93, [description="Parameter for quadrature axis mutual inductance saturation modelling"]

        # Transformer internal parameters
        RTfoPu = RTfPu * (UNomHV / UBaseHV) ^ 2 * (SNom / SnTfo), [description="Resistance of the generator transformer in pu (base SNom, UNom)"]
        XTfoPu = XTfPu * (UNomHV / UBaseHV) ^ 2 * (SNom / SnTfo), [description="Reactance of the generator transformer in pu (base SNom, UNom)"]
        # HACK: hardcoed one case
        # rTfoPu = ifelse(RTfPu > 0.0 || XTfPu > 0.0, UNomHV / UBaseHV / (UNomLV / UBaseLV), 1.0), [description="Ratio of the generator transformer in pu (base UBaseHV, UBaseLV)"]
        rTfoPu = UNomHV / UBaseHV / (UNomLV / UBaseLV), [description="Ratio of the generator transformer in pu (base UBaseHV, UBaseLV)"]
    end
end
@mtkmodel DynawoMachine begin
    @extend DynawoMachineParameterBase()
    @components begin
        terminal = Terminal()
        systembase = SystemBase()
        # ωRefPu = RealInput(), [input=true, description="Reference frequency in pu"]
        # PmPu = RealInput(), [input=true, description="Mechanical power in pu (base PNomTurb)"]
        # efdPu = RealInput(), [input=true, description="Input voltage of exciter winding in pu (user-selected base voltage)"]
        ωRefPu = RealInput()
        PmPu = RealInput()
        efdPu = RealInput()
    end
    @parameters begin
        # internal parameters for generator
        # Notation: Ra (resistance) + P ("'" or "Prim") + Pu (Per unit)
        Ra′Pu=0.003, [description="Armature resistance in pu"]
        Ld′Pu=0.15, [description="Direct axis stator leakage in pu"]
        Md′Pu=1.66, [description="Direct axis mutual inductance in pu"]
        LD′Pu=0.16634, [description="Direct axis damper leakage in pu"]
        RD′Pu=0.03339, [description="Direct axis damper resistance in pu"]
        Mrc′Pu=0, [description="Canay's mutual inductance in pu"]
        Lf′Pu=0.1699, [description="Excitation winding leakage in pu"]
        Rf′Pu=0.00074, [description="Excitation winding resistance in pu"]
        Lq′Pu=0.15, [description="Quadrature axis stator leakage in pu"]
        Mq′Pu=1.61, [description="Quadrature axis mutual inductance in pu"]
        LQ1′Pu=0.92815, [description="Quadrature axis 1st damper leakage in pu"]
        RQ1′Pu=0.00924, [description="Quadrature axis 1st damper resistance in pu"]
        LQ2′Pu=0.12046, [description="Quadrature axis 2nd damper leakage in pu"]
        RQ2′Pu=0.02821, [description="Quadrature axis 2nd damper resistance in pu"]
        MsalPu=0.05, [description="Constant difference between direct and quadrature axis saturated mutual inductances in pu"]
        # pu factor for excitation voltage
        Md′PuEfd, [description="Direct axis mutual inductance used to determine the excitation voltage in pu"]
        Md′PuEfdNom, [description="Direct axis mutual inductance used to determine the excitation voltage in nominal conditions in pu"]
    end
    @variables begin
        # # Input variables
        # ωRefPu(t), [input=true, description="Reference frequency in pu"]
        # PmPu(t), [input=true, description="Mechanical power in pu (base PNomTurb)"]
        # efdPu(t), [input=true, description="Input voltage of exciter winding in pu (user-selected base voltage)"]

        # Output variables
        ωPu(t), [output=true, description="Angular frequency in pu"]

        # d-q axis pu variables (base UNom, SNom)
        udPu(t), [description="Voltage of direct axis in pu"]
        uqPu(t), [description="Voltage of quadrature axis in pu"]
        idPu(t), [description="Current of direct axis in pu"]
        iqPu(t), [description="Current of quadrature axis in pu"]
        iDPu(t), [description="Current of direct axis damper in pu"]
        iQ1Pu(t), [description="Current of quadrature axis first damper sn pu"]
        iQ2Pu(t), [description="Current of quadrature axis second damper in pu"]
        ifPu(t), [description="Current of excitation winding in pu"]
        ufPu(t), [description="Voltage of exciter winding in pu (base voltage as per Kundur)"]
        λ_dPu(t), [description="Flux of direct axis in pu"]
        λ_qPu(t), [description="Flux of quadrature axis in pu"]
        λ_DPu(t), [description="Flux of direct axis damper in pu"]
        λ_fPu(t), [description="Flux of excitation winding in pu"]
        λ_Q1Pu(t), [description="Flux of quadrature axis 1st damper in pu"]
        λ_Q2Pu(t), [description="Flux of quadrature axis 2nd damper in pu"]

        # Other variables
        θ(t), [description="Rotor angle: angle between machine rotor frame and port phasor frame"]
        cmPu(t), [description="Mechanical torque in pu (base PNomTurb/ωNom)"]
        cePu(t), [description="Electrical torque in pu (base SNom/ωNom)"]
        PePu(t), [description="Electrical active power in pu (base SNom)"]

        # Saturated mutual inductances and related variables
        MdSat′Pu(t), [description="Direct axis saturated mutual inductance in pu"]
        MqSat′Pu(t), [description="Quadrature axis saturated mutual inductance in pu"]
        λ_AirGa′Pu(t), [description="Total air gap flux in pu"]
        λ_ADPu(t), [description="Common flux of direct axis in pu"]
        λ_AQPu(t), [description="Common flux of quadrature axis in pu"]
        mdsPu(t), [description="Direct axis saturated mutual inductance in the case when the total air gap flux is aligned on the direct axis in pu"]
        mqsPu(t), [description="Quadrature axis saturated mutual inductance in the case when the total air gap flux is aligned on the quadrature axis in pu"]
        cos2Eta(t), [description="Common flux of direct axis contribution to the total air gap flux in pu"]
        sin2Eta(t), [description="Common flux of quadrature axis contribution to the total air gap flux in pu"]
        miPu(t), [description="Intermediate axis saturated mutual inductance in pu"]
    end
    @equations begin
        # Park's transformations
        terminal.u_r ~ sin(θ) * udPu + cos(θ) * uqPu;
        terminal.u_i ~ (-cos(θ) * udPu) + sin(θ) * uqPu;
        terminal.i_r * systembase.SnRef / SNom ~ sin(θ) * idPu + cos(θ) * iqPu;
        terminal.i_i * systembase.SnRef / SNom ~ (-cos(θ) * idPu) + sin(θ) * iqPu;
        # Flux linkages
        λ_dPu ~ (MdSat′Pu + Ld′Pu + XTfoPu) * idPu + MdSat′Pu * ifPu + MdSat′Pu * iDPu;
        λ_fPu ~ MdSat′Pu * idPu + (MdSat′Pu + Lf′Pu + Mrc′Pu) * ifPu + (MdSat′Pu + Mrc′Pu) * iDPu;
        λ_DPu ~ MdSat′Pu * idPu + (MdSat′Pu + Mrc′Pu) * ifPu + (MdSat′Pu + LD′Pu + Mrc′Pu) * iDPu;
        λ_qPu ~ (MqSat′Pu + Lq′Pu + XTfoPu) * iqPu + MqSat′Pu * iQ1Pu + MqSat′Pu * iQ2Pu;
        λ_Q1Pu ~ MqSat′Pu * iqPu + (MqSat′Pu + LQ1′Pu) * iQ1Pu + MqSat′Pu * iQ2Pu;
        λ_Q2Pu ~ MqSat′Pu * iqPu + MqSat′Pu * iQ1Pu + (MqSat′Pu + LQ2′Pu) * iQ2Pu;
        # Equivalent circuit equations in Park's coordinates
        udPu ~ (Ra′Pu + RTfoPu) * idPu - ωPu * λ_qPu;
        uqPu ~ (Ra′Pu + RTfoPu) * iqPu + ωPu * λ_dPu;
        ufPu ~ Rf′Pu * ifPu + Dt(λ_fPu) / systembase.ωNom;
        0 ~ RD′Pu * iDPu + Dt(λ_DPu) / systembase.ωNom;
        0 ~ RQ1′Pu * iQ1Pu + Dt(λ_Q1Pu) / systembase.ωNom;
        0 ~ RQ2′Pu * iQ2Pu + Dt(λ_Q2Pu) / systembase.ωNom;
        # Mechanical equations
        Dt(θ) ~ (ωPu - ωRefPu.u) * systembase.ωNom;
        2 * H * Dt(ωPu) ~ cmPu * PNomTurb / SNom - cePu - DPu * (ωPu - ωRefPu.u);
        cePu ~ λ_qPu * idPu - λ_dPu * iqPu;
        PePu ~ cePu * ωPu;
        PmPu.u ~ cmPu * ωPu;
        # Excitation voltage pu conversion
        # ufPu ~ efdPu.u * (Kuf * rTfoPu);
        # HACK: fixed to excitation type "noload"
        #=
        Kuf=if ExcitationPu == ExcitationPuType.Kundur
            1
        elseif ExcitationPu == ExcitationPuType.UserBase
            Rf′Pu / Md′PuEfd
        elseif ExcitationPu == ExcitationPuType.NoLoad
            Rf′Pu / Md′Pu
        elseif ExcitationPu == ExcitationPuType.NoLoadSaturated
            Rf′Pu * (1 + md) / Md′Pu
        else
            Rf′Pu / Md′PuEfdNom
        end;
        =#
        ufPu ~ efdPu.u * (Rf′Pu / Md′Pu * rTfoPu);
        # Mutual inductances saturation
        λ_ADPu ~ MdSat′Pu * (idPu + ifPu + iDPu);
        λ_AQPu ~ MqSat′Pu * (iqPu + iQ1Pu + iQ2Pu);
        λ_AirGa′Pu ~ sqrt(λ_ADPu ^ 2 + λ_AQPu ^ 2);
        mdsPu ~ Md′Pu / (1 + md * λ_AirGa′Pu ^ nd);
        mqsPu ~ Mq′Pu / (1 + mq * λ_AirGa′Pu ^ nq);
        cos2Eta ~ λ_ADPu ^ 2 / λ_AirGa′Pu ^ 2;
        sin2Eta ~ λ_AQPu ^ 2 / λ_AirGa′Pu ^ 2;
        miPu ~ mdsPu * cos2Eta + mqsPu * sin2Eta;
        MdSat′Pu ~ miPu + MsalPu * sin2Eta;
        MqSat′Pu ~ miPu - MsalPu * cos2Eta;
    end
end

export Bus
# @mtkmodel Bus begin
#     @structural_parameters begin
#         N
#     end
#     @components begin
#         terminal = [Terminal() for i in 1:N]
#     end
#     @variables begin
#         u_r(t), [description="bus d-voltage", output=true]
#         u_i(t), [description="bus q-voltage", output=true]
#         i_r(t), [description="bus d-current", input=true]
#         i_i(t), [description="bus d-current", input=true]
#     end
#     @equations begin
#         [u_r ~ terminal[i].u_r for i in 1:N]...
#         [u_i ~ terminal[i].u_i for i in 1:N]...
#         i_r ~ sum([terminal[i].i_r for i in 1:N])
#         i_i ~ sum([terminal[i].i_i for i in 1:N])
#     end
# end

# @connector BusTerminal begin
#     u_r(t), [description="d-voltage", output=true]
#     u_i(t), [description="q-voltage", output=true]
#     i_r(t), [description="d-current", connect=Flow, input=true]
#     i_i(t), [description="q-current", connect=Flow, input=true]
# end

function Bus(_injectors...; name=:bus)
    injectors = Tuple[]
    for injector in _injectors
        if injector isa Tuple
            push!(injectors, injector)
        else
            push!(injectors, (injector, injector.terminal))
        end
    end

    # N = length(injectors)

    @variables begin
        u_r(t), [description="bus d-voltage", output=true]
        u_i(t), [description="bus q-voltage", output=true]
        i_r(t), [description="bus d-current", input=true]
        i_i(t), [description="bus d-current", input=true]
    end
    eqs = [
        [u_r ~ inj[2].u_r for inj in injectors]...
        [u_i ~ inj[2].u_i for inj in injectors]...
        i_r ~ sum([inj[2].i_r for inj in injectors])
        i_i ~ sum([inj[2].i_i for inj in injectors])
    ]
    ODESystem(eqs, t; systems=first.(injectors), name)
    # @named terminal = BusTerminal()
    # eqs = [connect(terminal, inj[2]) for inj in injectors]
    # ODESystem(eqs, t; systems=[terminal, first.(injectors)...], name)
end



# @mtkmodel DQBus begin
#     @variables begin
#         u_r(t), [description="d-voltage", output=true]
#         u_i(t), [description="q-voltage", output=true]
#         i_r(t), [description="d-current", input=true]
#         i_i(t), [description="d-current", input=true]
#     end
# end

# rotm(θ)=[cos(θ) -sin(θ); sin(θ) cos(θ)]
# @mtkmodel DQSwing begin
#     @extend DQBus()
#     @variables begin
#         ω(t)=0.0, [description="Rotor frequency"]
#         θ(t)=0.0, [description="Rotor angle"]
#         Pel(t), [description="Electrical Power"]
#     end
#     @parameters begin
#         M=1, [description="Inertia"]
#         D=0.1, [description="Damping"]
#         Pmech, [description="Mechanical Power"]
#         V=1.0, [description="Voltage magnitude"]
#     end
#     @equations begin
#         Dt(θ) ~ ω
#         Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
#         Pel ~ real(Complex(u_r, u_i) * conj(Complex(i_r, i_i)))
#         [u_r, u_i] ~ rotm(θ) * [V; 0]
#     end
# end

# @mtkmodel DQLine begin
#     @variables begin
#         src_u_r(t), [description="src d-voltage", input=false]
#         src_u_i(t), [description="src q-voltage", input=false]
#         dst_u_r(t), [description="dst d-voltage", input=false]
#         dst_u_i(t), [description="dst q-voltage", input=false]
#         src_i_r(t), [description="src d-current", output=true]
#         src_i_i(t), [description="src d-current", output=true]
#         dst_i_r(t), [description="dst d-current", output=true]
#         dst_i_i(t), [description="dst d-current", output=true]
#     end
# end

# @mtkmodel DQPiLine begin
#     @extend DQLine()
#     @parameters begin
#         R=1.0, [description="Resistance"]
#         X=1.0, [description="Reactance"]
#         src_B=0.0, [description="Shunt susceptance at src end"]
#         dst_B=0.0, [description="Shunt susceptance at dst end"]
#     end
#     @equations begin
#         src_i_r ~ -real(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
#         src_i_i ~ -imag(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
#         dst_i_r ~ -real(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
#         dst_i_i ~ -imag(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
#     end
# end