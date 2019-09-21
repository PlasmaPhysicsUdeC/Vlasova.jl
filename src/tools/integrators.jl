export verlet_position, verlet_velocity, mclachlan_velocity, BABAB, ChinA, ChinC, ABABABA, ABACABA

"""
Verlet scheme starting with position.
"""
verlet_position = VlasovaIntegrator("ABA", [0.5,
                                            1.0,
                                            0.5])

"""
Verlet scheme starting with velocity.
"""
verlet_velocity = VlasovaIntegrator("BAB", [0.5,
                                            1.0,
                                            0.5])

"""
Doc this! [`McLachlan et al. (1992)`](https://iopscience.iop.org/article/10.1088/0951-7715/5/2/011)
"""
mclachlan_velocity = VlasovaIntegrator("BABA", [1 - 1/sqrt(2),
                                                1/sqrt(2),
                                                1/sqrt(2),
                                                1 - 1/sqrt(2)])
"""
"""
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
                          gradient_coefficients = [ 1/48 ] )

ChinC = VlasovaIntegrator("ABACABA", [1/6,
                                      3/8,
                                      1/3,
                                      1/4,
                                      1/3,
                                      3/8,
                                      1/6],
                          gradient_coefficients = [ 1 / 48 ])

let
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


let
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
    # TODO: Check gradient coef
end



# TODO: Implement the **120 Eff gradient integrator from Omelyan2003
