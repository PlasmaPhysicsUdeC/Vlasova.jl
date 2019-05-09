export verlet_position, verlet_velocity, mclachlan_velocity

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
