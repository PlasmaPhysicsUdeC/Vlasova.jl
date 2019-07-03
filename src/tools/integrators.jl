export verlet_position, verlet_velocity, mclachlan_velocity, BABAB

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

# TODO: Implement the **120 Eff integrator from Omelyan2003
