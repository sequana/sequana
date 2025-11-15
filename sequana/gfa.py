class GFA:
    def __init__(self, filename):
        self.filename = filename
        self.read()

    def read(self):
        with open(self.filename, "r") as fin:
            for line in fin.readlines():
                if line.startswith("S"):
                    items = line.split()
                    assert items[0] == "S"
                    name = items[1]
                    sequence = items[2]
                    length = len(sequence)
                    depth = None
                    for item in items[2]:
                        if item.startswith("dp:i"):
                            depth = item.split(":")[-1]
                    print(name, length, depth)

                elif line.startswith("L"):
                    items = line.split()
                    assert items[0] == "L"
                    name = items[1]
                    edges = items[2]

                    for edge in edges.split(","):
                        ename = edge[0:-1]
                        direction = edge[-1]
                        print(f"L\t{name} {direction} {ename} 0M")
                else:
                    print(line)


# S edge11 ACGT
# L edge1 +  edge2 - 0M TC:i:15
#
# P contig1 edge1+,ede2- *
