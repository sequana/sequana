"""

We will not use SAM format since it is uncompressed. However,
a few tools may be added. Once could ba to convert SAM to BAM, 
which may already be present in other tools such as as pysam.


"""

import pandas as pd
import pylab
import pysam

# pysam uses htlib behing the scene and is very fast
# pysam works great for BAM file but with SAM, it needs to read the file after
# each compete iteration, which is not very useful


class SAM(pysam.AlignmentFile):
    """

    Header of a SAM file has N lines starting with '@' character.
    Using comment='@' and header None does not make the job.Using skiprows=N
    works but requires to know the number of lines in the header.


    .. todo:: read by chunk size for large files ?
    .. todo:: read two files for comparison ?


    FLAGS SAM format::

        1        0x1     template having multiple segments in sequencing
        2        0x2     each segment properly aligned according to the aligner
        4        0x4     segment unmapped
        8        0x8     next segment in the template unmapped
        16      0x10     SEQ being reverse complemented
        32      0x20     SEQ of the next segment in the template being reverse complemented
        64      0x40     the first segment in the template
        128      0x80    the last segment in the template
        256     0x100    secondary alignment
        512     0x200    not passing filters, such as platform/vendor quality controls
        1024     0x400   PCR or optical duplicate
        2048     0x800   supplementary alignme



    """
    def __init__(self, filename, *args):
        # works for py27 but not py35 probably a missing __init__  or __new__ in
        # AlignmentFile class. See e.g., stackoverflow/questions/
        # 26653401/typeerror-object-takes-no-parameters-but-only-in-python-3
        #super(SAM, self).__init__(filename, mode, *args)
        pysam.AlignmentFile.__init__(filename, *args)
        self.skiprows = self._guess_header_length()
        print('deprecated')

    def _guess_header_length(self):
        with open(self.filename, 'r') as fin:
            skiprows = 0
            while True:
                line = fin.readline()
                if line.startswith('@'):
                    skiprows += 1
                else:
                    break
        return skiprows

    def get_read_names(self):
        return self._get_column(0)

    def _get_column(self, col):
        data = []
        with open(self.filename, "r") as fin:
            for header in range(self.skiprows):
                fin.readline()

            for line in fin:
                data.append(line.split()[col])
        return data

    def plot_mapq_distribution(self):
        """Plot distribution of MAPQ scores (fifth column of SAM)


        The maximum MAPQ value that Bowtie 2 generates is 42 (though it doesn't
        say this anywhere in the documentation). In contrast, the maximum MAPQ
        value that BWA will generate is 37 (though once again, you -
        frustratingly - won't find this information in the manual).

        :reference: http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
        """
        pylab.clf()
        data = [float(x) for x in self._get_column(3)]
        pylab.hist(data, bins=100, normed=True)
        pylab.grid(True)
        pylab.xlabel('MAPQ score')
        pylab.ylabel('Fraction of reads')




