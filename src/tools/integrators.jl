export verlet_position, verlet_velocity, mclachlan_velocity, BABAB, ABABABA

verlet_position = VlasovaIntegrator("ABA", [0.5,
                                            1.0,
                                            0.5])
verlet_velocity = VlasovaIntegrator("BAB", [0.5,
                                            1.0,
                                            0.5])

mclachlan_velocity = VlasovaIntegrator("BABA", [1 - 1/sqrt(2),
                                                1/sqrt(2),
                                                1/sqrt(2),
                                                1 - 1/sqrt(2)])

BABAB = VlasovaIntegrator("BABAB", [0.19318332750378361,
                                    0.5,
                                    0.6136333449924327,
                                    0.5,
                                    0.19318332750378361 ])

begin
    lambda = 1/(2 - 2^(1/3))
    teta = lambda / 2

    global ABABABA = VlasovaIntegrator("ABABABA", [teta,
                                                   lambda,
                                                   0.5 - teta,
                                                   1 - 2*lambda,
                                                   0.5 - teta,
                                                   lambda,
                                                   teta])
end

# TODO: Implement the **120 Eff integrator from Omelyan2003
