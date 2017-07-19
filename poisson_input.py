def poisson_input(active, N_input, r, time_input, t1, t2):
    "function to built a poisson input"
    "with no zero firing rates"
    from brian import PoissonGroup,SpikeMonitor,reinit,clear,run,Hz,second
    import random as pyrandom

    reinit(states = True)
    clear(erase   = True, all = True)

    N  = 1000
    P  = PoissonGroup(N)
    S = SpikeMonitor(P)

    run(t1)

    P.rate = r
    run(time_input)

    P.rate = 0*Hz
    run(t2)

    s = []
    remove = []
    for i in xrange(N):
        s.append(len(S[i]))
        if len(S[i]) == 0:
            remove.append(i)

    pos = [x for x in range(N) if x not in remove]
    pos = pyrandom.sample(pos, len(active))

    C =  []
    for ii in xrange(len(pos)):
        D = list(S.spiketimes[pos[ii]])
        D = [(active[ii],x) for x in D]
        C += D

    C = [(i,t *second) for (i,t) in C]
    C = sorted(C, key=lambda x: x[1])

    reinit(states = True)
    clear(erase   = True, all = True)

    return C
