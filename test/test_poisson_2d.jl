using Test, Vlasova
import FFTW

# Make plasma
box = Box( name = "test",
           Nx = (8, 10),
           Nv = (12, 14),
           Lx = (5pi, 10.0),
           vmin = (-6, -8),
           vmax = (5, 7)
           )

electrons = Specie(name = "electrons",
                   charge = -1.0,
                   mass = 1.0,
                   temperature = 1.0,
                   distribution = ones(box.Nx...) ⊗ maxwellian2d( box.v... )  )

plasma = Plasma(electrons, box);

# Time info
dt = 0.1
final_time = 100

# Make structs
integrator = verlet_velocity

poisson = Vlasova.Poisson( plasma )

sadv = Vlasova.SpaceAdvection( plasma,
                       integrator,
                       dt )

vadv = Vlasova.VelocityAdvection( plasma,
                                  integrator,
                                  dt,
                                  FFTW_flags = FFTW.ESTIMATE)

tm = Vlasova.TimeManager( final_time,
                          dt )

# Tests structs
Nx2 = div(box.Nx[1], 2)
Ny2 = div(box.Nx[2], 2)

kx = Array(0:Nx2) * 2pi / box.Lx[1]
ky = FFTW.fftshift( Array(1:box.Nx[2]) .- (Ny2 + 1) ) * 2pi / box.Lx[2]
k2 = kx.^2 .+ ky'.^2

## test Poisson
@test poisson.k[1] == (kx .* ones(box.Nx[2])')
@test poisson.k[2] == (ones( Nx2 + 1 ) .* ky')
@test poisson.pot2dens == k2

let one_over_k2 = 1 ./ k2
    one_over_k2[1] = 0.0
    @test poisson.dens2field[1] == ( -1im * kx  .* one_over_k2 )
    @test poisson.dens2field[2] == ( -1im * ky' .* one_over_k2 )
end

nothing
