module BEMfunctions
      import LinearAlgebra; la = LinearAlgebra
      using rotationalMatrices
      using BEMclasses

      export positions, angles, forces, inflowQS, inflowDyn

      ρ = 1.225

      function Cl(al)
            rad2deg(al)
            return 0.22*al
      end

      function Cd(al)
            return 0.25
      end

      function fg(a)
            if a <= 0.33
                  return 1
            else
                  return 0.25*(5-3a)
            end
      end




      function positions(blade :: blade, WTG :: turbine)
            a4b = zRot(blade.θr)
            rt = [WTG.H0;0.;0.]
            rs = [0.;0.;-WTG.Ls]

            blade.pos = rt.+WTG.a01*WTG.a12*rs.+
                        WTG.a01*WTG.a12*WTG.a34*a4b*blade.r
            return nothing
      end

      function angles(blade :: blade)
            blade.ϕ[:] = atan.(blade.Vb[3,:], -blade.Vb[2,:])
            blade.α[:] = blade.ϕ[:]-(blade.β[:] .+ blade.θp)
            return nothing
      end

      function forces(blade :: blade)
            V2 = sum(blade.Vb.^2,dims=1)

            blade.L[:] = 0.5*ρ*V2.*blade.c.* Cl.(blade.α)
            blade.D[:] = 0.5*ρ*V2.*blade.c.* Cd.(blade.α)
            return nothing
      end

      function inflowQS(blade :: blade, WTG :: turbine)

            a = -blade.w[3:3,:] ./ blade.V0[3:3,:]
            # Prandtl's tip loss correction:
            f = WTG.B*0.5 * (blade.R.-blade.r[1:1,:])./
                (blade.r[1:1,:].*abs.(sin.(blade.ϕ)))
            F = 2/pi*acos.(exp.(-f))

            gVel = sqrt.(blade.V0[2:2,:].^2+
                         (blade.V0[3:3,:]+fg.(a).*blade.wqs[3:3,:]).^2)

            blade.wqs[2,:] = -WTG.B*blade.L.*sin.(blade.ϕ)./
                              (4pi*ρ*blade.r[1:1,:]./blade.R.*F.*gVel)

            blade.wqs[3,:] = -WTG.B*blade.L.*cos.(blade.ϕ)./
                              (4pi*ρ*blade.r[1:1,:]./blade.R.*F.*gVel)


            return nothing
      end

      function inflowDyn(blade :: blade, WTG :: turbine, Δt)
            a = -blade.w[3:3,:] ./ blade.V0[3:3,:]
            # Prandtl's tip loss correction:
            f = WTG.B*0.5 * (blade.R.-blade.r[1:1,:])./
                (blade.r[1:1,:].*abs.(sin.(blade.ϕ)))
            F = 2/pi*acos.(exp.(-f))

            gVel = sqrt.(blade.V0[2:2,:].^2+
                         (blade.V0[3:3,:]+fg.(a).*blade.wqs[3:3,:]).^2)

            wqsOld = blade.wqs

            blade.wqs[2,:] = -WTG.B*blade.L.*sin.(blade.ϕ)./
                              (4pi*ρ*blade.r[1:1,:]./blade.R.*F.*gVel)

            blade.wqs[3,:] = -WTG.B*blade.L.*cos.(blade.ϕ)./
                              (4pi*ρ*blade.r[1:1,:]./blade.R.*F.*gVel)


            τ1 = (1.1*blade.R)./(blade.V0[3:3,:].+1.3*blade.w[3:3,:])
            τ2 = τ1.*(0.39.-0.26*(blade.r[1:1,:]./blade.R).^2)

            H = blade.wqs + 0.6*τ1.*(blade.wqs.-wqsOld)/Δt

            blade.wint = H .+ (blade.wint.-H).*exp.(-Δt./τ1)

            blade.w = blade.wint .+ (blade.w-blade.wint).*exp.(-Δt./τ2)
            return nothing
      end


end
