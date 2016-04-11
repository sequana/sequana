"""

USed for adapter removal benchmarks. 

May not be used on long term. 

May be removed

"""


class AdapterRemoval(object):
    def __init__(self):
        pass


def adapters_files_to_list(filename1, filename2):

    fh1 = open(filename1, 'r')
    fh2 = open(filename1, 'r')
    data1 = fh1.readlines()
    data2 = fh2.readlines()
    fh1.close()
    fh2.close()

    len(data1) == len(data2), "incompatible files. Must have same length"

    fh = open("adapters_list.fa", 'w')
    count = 0
    for line1, line2 in zip(data1, data2):
        line1 = line1.strip()
        line2 = line2.strip()
        if line1.startswith(">"):
            pass
        else:
            fh.write(line1+" " +line2+ "\n")
            count += 1

    fh.close()
    print("Saved %s adapters in adapters_combined.fa" % count)


def adapters_to_clean_ngs(filename):
    fh1 = open(filename, 'r')
    data1 = fh1.readlines()
    fh1.close()

    count = 0
    fh = open("adapters_ngs.txt", "w")
    for line in data1:
        line = line.strip().strip("\n")
        if line.startswith('>'):
            pass
        else:
            data = "adapter_%s\t%s\t0.5\t31\t10\t0\t0\n"% (count+1, line)
            fh.write(data)
            count+=1
    fh.close()



adapters_to_clean_ngs("adapters_48_PCR-free_FWD.fa")
#if __name__ == "__main__":
#    import sys
#    args = sys.argv
#    adapters_files_to_list(args[1], args[2])
