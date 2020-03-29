module BEMclasses

    export turbine, blade

    mutable struct turbine
        H0    :: Float64
        Ls    :: Float64
        R     :: Float64
        tR    :: Float64


        θtilt :: Float64
        θcone :: Float64
        θyaw  :: Float64
        a01   :: Array{Float64,2}
        a10   :: Array{Float64,2}
        a12   :: Array{Float64,2}
        a21   :: Array{Float64,2}
        a34   :: Array{Float64,2}
        a43   :: Array{Float64,2}

        name  :: String
        B     :: Int
        Thrust :: Float64
        Torque :: Float64
        Power  :: Float64
    end

    mutable struct blade
        nE :: Int
        R  :: Float64
        r  :: Array{Float64}
        c  :: Array{Float64}
        θr :: Float64
        θp :: Float64
        β  :: Array{Float64}
        α  :: Array{Float64}
        ϕ  :: Array{Float64}

        L  :: Array{Float64}
        D  :: Array{Float64}

        pos :: Array{Float64}
        V0 :: Array{Float64}
        Vb :: Array{Float64}

        # Induction parameters
        wqs  :: Array{Float64}
        wint :: Array{Float64}
        w    :: Array{Float64}

    end

end
