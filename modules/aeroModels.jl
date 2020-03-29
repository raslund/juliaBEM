module aeroModels

      import LinearAlgebra; la = LinearAlgebra

      export towerShadow, windShear

      function towerShadow(pos,U, H0, a=7.)
            if pos[1] > H0
                  return U
            else
                  Upot = la.norm(U)
                  r = la.norm(pos[2:3])
                  Ur = Upot*(1.0-(a/r)^2.0)*pos[3]/r
                  Uθ = -Upot*(1.0+(a/r)^2.0)*(-pos[2]/r)

                  U[3] = Ur*pos[3]/r - Uθ*(-pos[2]/r)
                  U[2] = -Ur*(-pos[2]/r) - Uθ*pos[3]/r

                  return U
            end
      end

      function windShear(pos,U, H0, nu=0.2)
            Uout = U.*(pos[1]/H0).^nu
            return Uout
      end
end
