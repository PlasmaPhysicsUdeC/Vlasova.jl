export verlet_position, verlet_velocity, mclachlan_velocity, BABAB, ChinA, ChinC, ABABABA, ABACABA

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

ChinA = VlasovaIntegrator("BACAB", [1/6,
                                    1/2,
                                    2/3,
                                    1/2,
                                    1/6],
                          gradient_coefficients = [2 / (3*48) ])

ChinC = VlasovaIntegrator("ABACABA", [1/6,
                                      3/8,
                                      1/3,
                                      1/4,
                                      1/3,
                                      3/8,
                                      1/6],
                          gradient_coefficients = [1 / (4*48) ])

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


begin
    xi = 0.006938106540706989
    lambda = 0.2470939580390842
    teta = 0.08935804763220157
    
    global ABACABA = VlasovaIntegrator("ABACABA", [teta,
                                                   lambda,
                                                   0.5 - teta,
                                                   1 - 2*lambda,
                                                   0.5 - teta,
                                                   lambda,
                                                   teta],
                                       gradient_coefficients = [xi])
end



# TODO: Implement the **120 Eff integrator from Omelyan2003
