import PyPlot; plt = PyPlot
import LinearAlgebra; la = LinearAlgebra
import Trapz
import Random
using ProgressMeter

# Add path to custom modules:
push!(LOAD_PATH, pwd()*"/modules/")
using rotationalMatrices
using aeroModels
using BEMclasses
using statusPrints
using BEMfunctions

startup()

# Turbine initialization
H0 = 119.0; Ls = 7.1; R = 89.15; towerRadius = 3.5
θtilt = deg2rad(-5.2); θcone = deg2rad(1.0);
θyaw = deg2rad(0.0)

# Rotational matrices
# Yaw matrix
a01 = xRot(θyaw)
a10 = la.transpose(a01)

# Tilt matrix
a12 = yRot(θtilt)
a21 = la.transpose(a12)

# Cone matrix
a34 = yRot(θcone)
a43 = la.transpose(a34)

WTG = turbine(H0,
              Ls,
              R,
              towerRadius,
              θtilt,
              θcone,
              θyaw,
              a01,
              a10,
              a12,
              a21,
              a34,
              a43,
              "DTU10MW",
              3,
              0.,
              0.,
              0.)


rt = [H0;0.;0.]
rs = [0.;0.;-Ls]

# Blade initialization
numElements = 15;
rb0 = Array{Float64}(undef,(3,numElements))
rb0[1,:] = range(WTG.R/numElements,stop=0.99*WTG.R,length=numElements)

θb0 = -pi/2.0
θp = 0
β = zeros(1,numElements)

# Initialize turbine blades
blades = Array{blade, 1}(undef, WTG.B)

for i in 1:length(blades)
      blades[i] = blade(numElements,
                        WTG.R,
                        rb0,
                        ones((1,numElements)),
                        -(θb0 + (i-1)*(2.0*pi/3)),
                        θp,
                        β,
                        zeros((1,numElements)),
                        zeros((1,numElements)),
                        zeros((1,numElements)),
                        zeros((1,numElements)),
                        rt .+ a01*a12*rs .+ a01*a12*a34*
                           zRot(-(θb0 + (i-1)*(2.0*pi/3)))*rb0,
                        (zeros((3,numElements))),
                        (zeros((3,numElements))),
                        (zeros((3,numElements))),
                        (zeros((3,numElements))),
                        (zeros((3,numElements)))
                        )
end

# Operational specifics:
U0 = [0;0;10.0]
ρ = 1.225
ω0 = -0.62

# Simulation time
Ttotal = 600
Δt = 0.02
Nt = ceil(Int,Ttotal/Δt)+1
Time = range(0,length=Nt,step=Δt)

# # Threads.@threads

function fg(a)
      if a <= 0.33
            return 1
      else
            return 0.25*(5-3a)
      end
end



pBar = Progress(Nt-1, dt = 0.2,
             barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=30)

if Threads.nthreads() > 1
      message(["Running on $(Threads.nthreads()) threads!"])
end


@timev for t in 1:Nt-1

      WTG.Thrust = 0
      WTG.Torque = 0

      Threads.@threads for b in 1:WTG.B
            blades[b].θr = blades[b].θr + ω0 * Δt

            # Blade element position in S0
            positions(blades[b],WTG)
            # blades[b].V0 .= 1.0*U0

            for e in 1:blades[b].nE
                  a4b = zRot(blades[b].θr)
                  ab4 = la.transpose(a4b)

                  # Wind velocity in S0 (Wind models)
                  V0 = 1.0*U0
                  V0 = windShear(blades[b].pos[:,e],U0, WTG.H0)
                  V0 = towerShadow(blades[b].pos[:,e],V0, WTG.H0, WTG.tR)
                  blades[b].V0[:,e] = V0
                  # Wind velocity in Sb (Velocity triangle)
                  Vb = ab4*a43*a21*a10*V0
                  Vrot = [0,ω0* (a34*blades[b].r[:,e])[1],0]
                  blades[b].Vb[:,e] = blades[b].w[:,e] + Vb + Vrot
            end

            # Angles and forces
            angles(blades[b])
            forces(blades[b])

            # Dynamic Inflow
            inflowDyn(blades[b],WTG,Δt)

      end

      # Integrated rotor quantities
      for b in 1:WTG.B
      WTG.Thrust += sum( blades[b].L[:].*cos.(blades[b].ϕ[:])+
                         blades[b].D[:].*sin.(blades[b].ϕ[:]))
      WTG.Torque += sum( blades[b].L[:].*sin.(blades[b].ϕ[:])-
                         blades[b].D[:].*cos.(blades[b].ϕ[:]))
      end
      # plt.figure("zpos theta")
      # plt.scatter(Time[t],WTG.Thrust, c="g")
      #

      next!(pBar)
end

# plt.display_figs()
