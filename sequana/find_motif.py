from sequana import BAM, FastQ



def find_motif(bamfile, motif, window=200):


    b1 = BAM(bamfile)

    # FIND motif and create pictures
    count = 0
    found = []
    for a in b1:
        count +=1
        if a.query_sequence is None:
            continue
        seq = a.query_sequence; 
        # Here we search for CAGCAGCAG no need to go further since there are errors
        # looking for CAG only leads to less robust threshold but could be tuned
        # further of course
        X1 = [seq[i:i+window].count(motif) for i in range(len(seq))]

        S = sum([x>5 for x in X1])

        if S > 10:
            print(a.query_name, S, a.reference_name)
            found.append(a)
            off = a.query_alignment_start
            clf()
            plot(range(off+a.reference_start, off+a.reference_start+len(seq)),X1)
            savefig("{}_{}_{}.png".format(a.reference_name, S, a.query_name.replace("/", "_")))

        if count%10000 == 0: 
            print(count)



