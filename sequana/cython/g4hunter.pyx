# cython: profile=True

cdef _cython_base_score(str line):
    cdef list[int] scores = []
    cdef int counterG = 0
    cdef int counterC = 0
    cdef int i
    cdef str item

    for item in line:
        if item in "Gg":
            if counterC:
                C = min(counterC, 4)
                for i in range(counterC):
                    scores.append(-1 * C)
            counterG += 1
            counterC = 0
        elif item in "Cc":
            if counterG:
                G = min(counterG, 4)
                for i in range(counterG):
                    scores.append(G)
            counterG = 0
            counterC += 1
        else:
            if counterG:
                G = min(counterG, 4)
                for i in range(counterG):
                    scores.append(G)
            if counterC:
                C = min(counterC, 4)
                for i in range(counterC):
                    scores.append(-1 * C)

            scores.append(0)
            counterC = 0
            counterG = 0
    G = min(counterG, 4)
    for i in range(counterG):
        scores.append(G)
    C = min(counterC, 4)
    for i in range(counterC):
        scores.append(-1 * C)
    return scores

def cython_base_score(line):
    return _cython_base_score(line)
