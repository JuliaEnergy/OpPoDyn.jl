@mtkmodel PiLine_fault begin
    @parameters begin
        R, [description="Resistance of branch in pu (base unclear?)"]
        X, [description="Reactance of branch in pu (base unclear?)"]
        G_src, [description="Conductance of src shunt (base unclear?)"]
        B_src, [description="Susceptance of src shunt (base unclear?)"]
        G_dst, [description="Conductance of dst shunt (base unclear?)"]
        B_dst, [description="Susceptance of dst shunt (base unclear?)"]
        r_src=1, [description="Src end transformation ratio"]
        r_dst=1, [description="Src end transformation ratio"]
        #R_fault=0, [description="Fault Resistance in pu (base unclear?)"]
        #X_fault=0, [description="Fault Reactance in pu (base unclear?)"]
        G_fault=0, [description="Fault Resistance in pu (base unclear?)"]
        B_fault=0, [description="Fault Reactance in pu (base unclear?)"]
        pos, [description="Fault Position (from src, percent of the line)"]
        active=1, [description="Line active or switched off"]
        shortcircuit=0, [description="shortcircuit on line"] #false, [description="shortcircuit on line"]
        #faultimp, [description="1 if fault impedance given, else 0"] #blöd, dass man dann auch extra angeben muss, dass keine Fehlerimpedanz da ist. Aber vertretbar. Wie kann man error ausgeben, wenn G_fault oder B_fault gegeben ist, aber faultimp=0 ist?
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    begin
        Z = R + im*X
        #Y_f = (G_fault + im * B_fault) * shortcircuit #1/(R_fault + im*X_fault)
        Z_a = Z * pos
        Z_b = Z * (1-pos)
        Ysrc = G_src + im*B_src
        Ydst = G_dst + im*B_dst
        Vsrc = src.u_r + im*src.u_i
        Vdst = dst.u_r + im*dst.u_i
        V₁ = r_src * Vsrc
        V₂ = r_dst * Vdst
        i₁ = Ysrc * V₁
        i₂ = Ydst * V₂
        V_mnormal = V₁*(1-pos) + V₂*pos #(V₁*Z_b + V₂*Z_a)/(Z_a+Z_b)
        V_m = V_mnormal * (1-shortcircuit) + shortcircuit * (0+0*im) #ifelse(shortcircuit>0, 0.0, V_mnormal)
        #V_m = (V₁*(1-pos) + V₂*pos)/(1 + Y_f * Z * pos * (1-pos)) * (shortcircuit * faultimp + (1 - shortcircuit)) #Problem: nicht das korrekte Ergebnis -> liegt das jetzt daran, dass keine Pi-Line verwendet wird? Jein, offenbar ändert das schon was, aber auch mit Pi Line ist es nicht richtig...
        V_a = V₁ - V_m
        V_b = V_m - V₂
        i_m2 = V_b / Z_b
        i_m1 = V_a / Z_a
        isrc = -i_m1 - i₁
        idst = i_m2 - i₂
        i_f = i_m1 - i_m2

        #Y_a = Ysrc * pos #Annahme: Ysrc=Ydst=Y_quer/2
        #Y_b = Ysrc * (1-pos)
        #i_11 = Y_a * V₁
        #i_22 = Y_b * V₂
        #V_mnormal = (V₁/Z_a + V₂/Z_b)/(Y_a+Y_b+1/Z_a+1/Z_b)
        #V_mnormal = (V₁ * (1-pos) + V₂ * pos)/(Ysrc * Z * pos * (1-pos) + 1) #Z und Y eingesetzt und zusammengefasst
        #V_m = V_mnormal * (1-shortcircuit) + shortcircuit * (0+0*im) #ifelse(shortcircuit>0, 0.0, V_mnormal)
        #i_12 = Y_a * V_m
        #i_21 = Y_b * V_m
        #V_a = V₁ - V_m
        #V_b = V_m - V₂
        #i_a = V_a / Z_a
        #i_b = V_b / Z_b
        #i_m1 = i_a - i_12
        #i_m2 = i_21 + i_b
        #isrc = -i_a - i_11
        #idst = i_b - i_22
        #i_f = i_m1 - i_m2
        

        #Z_ges = Z_a+ (Ysrc * Z * (1-pos)+1)/(Y_f + Ysrc * (2 + (1-pos) * Z * (Ysrc + Y_f))) #gesamter Widerstand über den V_1 abfällt (ohne Y_src), im pi Modell, vereinfacht
        #i_a = V₁ / Z_ges
        #V_a = i_a * Z_a
        #V_m = V₁ - V_a
        #i_f = V_m * Y_f
        #i_12 = V_m * Y_a
        #i_21 = V_m * Y_b
        #i_m1 = i_a - i_12
        #i_m2 = i_m1 - i_f
        #i_b = i_m2 - i_21
        #isrc = -i_a - i_11
        #idst = i_b - i_22


        #Z_ges = Z_a+ 1/(Y_f + 1/(Z_b+1/Ydst)) #gesamter Widerstand über den V_1 abfällt (ohne Y_src)
        #i_m1 = V₁ * 1/Z_ges
        #V_a = i_m1 * Z_a
        #V_fault = V₁ - V_a
        #V_b = V_fault - V₂
        #i_f = Y_f * V_fault
        #i_m2 = i_m1 - i_f
        #isrc = -i_m1 - i₁
        #idst =  i_m2 - i₂
    end
    @equations begin
        src.i_r ~ active * simplify(real(isrc))
        src.i_i ~ active * simplify(imag(isrc))
        dst.i_r ~ active * simplify(real(idst))
        dst.i_i ~ active * simplify(imag(idst))
    end
end
