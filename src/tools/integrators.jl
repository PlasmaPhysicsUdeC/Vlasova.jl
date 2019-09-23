export verlet_position,
    verlet_velocity,
    mclachlan_velocity,
    BABAB,
    chin_A,
    chin_C,
    ABABABA,
    ABACABA

"""
```julia
VlasovaIntegrator("ABA", [0.5,
                          1.0,
                          0.5])
```

Developed by [`Verlet (1967)`](https://journals.aps.org/pr/abstract/10.1103/PhysRev.159.98).
"""
verlet_position = VlasovaIntegrator("ABA", [0.5,
                                            1.0,
                                            0.5])

"""
```julia
VlasovaIntegrator("BAB", [0.5,
                          1.0,
                          0.5])
```

Developed by [`Swope et al. (1982)`](https://aip.scitation.org/doi/10.1063/1.442716).
"""
verlet_velocity = VlasovaIntegrator("BAB", [0.5,
                                            1.0,
                                            0.5])

"""
```julia
VlasovaIntegrator("BABA", [1 - 1/sqrt(2),
                           1/sqrt(2),
                           1/sqrt(2),
                           1 - 1/sqrt(2)])
```

Developed by [`McLachlan et al. (1992)`](https://iopscience.iop.org/article/10.1088/0951-7715/5/2/011).
"""
mclachlan_velocity = VlasovaIntegrator("BABA", [1 - 1/sqrt(2),
                                                1/sqrt(2),
                                                1/sqrt(2),
                                                1 - 1/sqrt(2)])
"""
```julia
VlasovaIntegrator("BABAB", [0.19318332750378361,
                            0.5,
                            0.6136333449924327,
                            0.5,
                            0.19318332750378361 ])
```
"""
BABAB = VlasovaIntegrator("BABAB", [0.19318332750378361,
                                    0.5,
                                    0.6136333449924327,
                                    0.5,
                                    0.19318332750378361 ])

"""
```julia
VlasovaIntegrator("BACAB", [ 1/6,
                             1/2,
                             2/3,
                             1/2,
                             1/6 ],
                          gradient_coefficients = [ 1/48 ] )
```

Developed by [`Chin (1997)`](https://www.sciencedirect.com/science/article/abs/pii/S0375960197000030).
"""
chin_A = VlasovaIntegrator("BACAB", [1/6,
                                     1/2,
                                     2/3,
                                     1/2,
                                     1/6],
                           gradient_coefficients = [ 1/48 ] )

"""
```julia
VlasovaIntegrator("ABACABA", [ 1/6,
                               3/8,
                               1/3,
                               1/4,
                               1/3,
                               3/8,
                               1/6 ],
                            gradient_coefficients = [ 1 / 48 ] )
```

Developed by [`Chin (1997)`](https://www.sciencedirect.com/science/article/abs/pii/S03759601970000308).
"""
chin_C = VlasovaIntegrator("ABACABA", [1/6,
                                       3/8,
                                       1/3,
                                       1/4,
                                       1/3,
                                       3/8,
                                       1/6],
                           gradient_coefficients = [ 1 / 48 ])

"""
```julia
VlasovaIntegrator("ABABABA", [ θ,
                               λ,
                               0.5 - θ,
                               1 - 2λ,
                               0.5 - θ,
                               λ,
                               θ])
```
where

```julia
λ = 1/(2 - 2^(1/3))
θ = λ / 2
```
"""
ABABABA = VlasovaIntegrator("ABABABA",
                                   begin
                                   λ = 1/(2 - 2^(1/3))
                                   θ = λ / 2
                                   [ θ,
                                     λ,
                                     0.5 - θ,
                                     1 - 2λ,
                                     0.5 - θ,
                                     λ,
                                     θ ]
                                   end)


"""
```julia
VlasovaIntegrator("ABACABA", [ θ,
                               λ,
                               0.5 - θ,
                               1 - 2λ,
                               0.5 - θ,
                               λ,
                               θ ],
                             gradient_coefficients = [ ξ ])
```

where

```julia
λ = 0.2470939580390842
θ = 0.08935804763220157
ξ = 0.006938106540706989
```
"""
ABACABA = VlasovaIntegrator("ABACABA",
                            begin
                            λ = 0.2470939580390842
                            θ = 0.08935804763220157
                            ξ = 0.006938106540706989
                            [ θ,
                              λ,
                              0.5 - θ,
                              1 - 2*λ,
                              0.5 - θ,
                              λ,
                              θ ]
                            end,
                            gradient_coefficients = [ ξ ])



# TODO: Implement the **120 Eff gradient integrator from Omelyan2003
