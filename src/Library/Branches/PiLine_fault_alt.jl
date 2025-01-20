@mtkmodel PiLine_fault begin
    @parameters begin
        R, [description="Resistance of branch in pu (base unclear?)"]
        X, [description="Reactance of branch in pu (base unclear?)"]
        R_fault=0, [description="Fault Resistance in pu (base unclear?)"]
        X_fault=0, [description="Fault Reactance in pu (base unclear?)"]
        pos, [description="Fault Position (from src, percent of the line)"]
        G_src, [description="Conductance of src shunt (base unclear?)"]
        B_src, [description="Susceptance of src shunt (base unclear?)"]
        G_dst, [description="Conductance of dst shunt (base unclear?)"]
        B_dst, [description="Susceptance of dst shunt (base unclear?)"]
        r_src=1, [description="Src end transformation ratio"]
        r_dst=1, [description="Src end transformation ratio"]
        active=1, [description="Line active or at fault"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    begin
        #R, X, R_fault, X_fault, pos, G_dst, B_dst = symbols("R X R_fault X_fault pos G_dst B_dst")
        #im = SymPy.I 
        Z = R + im*X
        Z_f = R_fault + im*X_fault
        Z_a = Z * pos
        Z_b = Z * (1-pos)
        Ysrc = G_src + im*B_src
        Ysrc_a = Ysrc * pos
        Ysrc_b = Ysrc * (1-pos)
        Ydst = G_dst + im*B_dst
        Ydst_a = Ydst * pos
        Ydst_b = Ydst * (1-pos)
        Vsrc = src.u_r + im*src.u_i
        Vdst = dst.u_r + im*dst.u_i
        V₁ = r_src * Vsrc
        V₂ = r_dst * Vdst
        i₁ = Ysrc * V₁
        i₂ = Ydst * V₂
        iₘ = 1/Z * (V₁ - V₂)
        isrc = -iₘ - i₁
        idst =  iₘ - i₂
    end
    @equations begin
        src.i_r ~ active * simplify(real(isrc))
        src.i_i ~ active * simplify(imag(isrc))
        dst.i_r ~ active * simplify(real(idst))
        dst.i_i ~ active * simplify(imag(idst))
    end
end
#function piline_fault(; R, X, B, pos, X_f=0, R_f=0, src, dst)
    #R_a = R * pos
    #X_a = X * pos
    #B_srca = pos * B/2
    #B_dsta = pos * B/2
    #@named pibranch_a = PiLine(;R=R_a, X=X_a, B_src=B_srca, B_dst=B_dsta, G_src=0, G_dst=0)#PiLine_fault(;R, X, B_src=B/2, B_dst=B/2, G_src=0, G_dst=0, pos)
    #R_b = R * (1-pos)
    #X_b = X * (1-pos)
    #B_srcb = (1-pos) * B/2
    #B_dstb = (1-pos) * B/2
    #@named pibranch_b = PiLine(;R=R_b, X=X_b, B_src=B_srcb, B_dst=B_dstb, G_src=0, G_dst=0)
    #@named mtkbusfault = LineFaultBus(; R_f=0, X_f=0)
   # @named bus_fault = Bus(mtkbusfault; vidx=0)
  #  @named l5f = Line(MTKLine(pibranch_a), src=src, dst=0)
 #   @named lf7 = Line(MTKLine(pibranch_b), src=0, dst=dst)
#end

#@named l57 = Line(piline_fault(; R=0.0320, X=0.1610, B=0.3060, pos=0.99), src=5, dst=7)
#piline_fault(; R=0.0320, X=0.1610, B=0.3060, pos=0.99, src=5, dst=7)