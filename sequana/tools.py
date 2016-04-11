import os
import string

from sequana import BAM


__all__ = ['bam_to_mapped_unmapped_fastq']


try:
    _translate = string.maketrans('ACGTacgt', 'TGCAtgca')
except:
    _translate = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')


def reverse_complement(seq):
    # use hastag but for now, do it manually
    return seq.translate(_translate)[::-1]


def bam_to_mapped_unmapped_fastq(filename, mode='pe'):
    """Create mapped and unmapped fastq files from a BAM file

    Given a BAM file, create FASTQ with R1/R2 reads mapped and unmapped.
    In the paired-end case, 4 files are created.

    The interest of this function is that it does not create intermediate files
    limiting IO in the process.

    Reads that are not paired or not "proper pair" (neither flag 4 nor flag 8)
    are ignored

    Secondary alignment (flag 256) are dropped so as to remove any ambiguous alignments and
    keep the number of final reads equal to the initial number of reads

    If R1 is mapped or R2 is mapped then the reads are considered mapped. If
    both R1 and R2 are unmapped, then reads are unmapped.

    Note about chimeric alginment: one is the representative and the other is
    the supplementary. This flag is not used in this function. Note also that
    chimeric alignment have same QNAME and flag 4 and 8

    .. todo:: comments is currently harcoded and should be removed. comments
        from the fastq qname are not stored in BAM.
    """
    bam = BAM(filename)

    newname, ext = os.path.splitext(filename)

    if mode == "pe":
        import collections
        stats = collections.defaultdict(int)
        unpaired = 0
        stats['R1_unmapped'] = 0
        stats['R2_unmapped'] = 0
        stats['R1_mapped'] = 0
        stats['R2_mapped'] = 0
        stats['duplicated'] = 0
        stats['unpaired'] = 0

        R1_mapped = open(newname + "_R1.mapped.fastq", "bw")
        R2_mapped = open(newname + "_R2.mapped.fastq", "bw")
        R1_unmapped = open(newname + "_R1.unmapped.fastq", "bw")
        R2_unmapped = open(newname + "_R2.unmapped.fastq", "bw")

        from easydev import Progress
        pb = Progress(len(bam))
        bam.reset()
        for i, this in enumerate(bam):
            # check that this is paired indeed
            if this.is_paired is False:
                stats['unpaired'] += 1
                # What to do in such case ?
            elif this.flag & 256:
                # Unmapped reads are in the BAM file but have no valid assigned
                # position (N.B., they may have an assigned position, but it should be ignored).
                # It's typically the case that a number of reads can't be aligned, due to things
                # like sequencing errors, imperfect matches between the DNA sequenced and the
                # reference, random e. coli or other contamination, etc..
                # A secondary alignment occurs when a given read could align reasonably well to
                # more than one place. One of the possible reported alignments is termed "primary"
                # and the others will be marked as "secondary".
                stats['secondary'] +=1 
            else:
                # inpysam, seq is a string and qual a bytes....
                if this.is_reverse is True:
                    #print("reversed complement entry %s" % i)
                    txt = b"@" + bytes(this.qname, "utf-8") + b"\t%s\n"
                    revcomp = reverse_complement(this.seq)
                    txt += bytes(revcomp, "utf-8") + b"\n"
                    txt += b"+\n"
                    txt += bytes(this.qual[::-1], 'utf-8') + b"\n"
                else:
                    txt = b"@" + bytes(this.qname, "utf-8") + b"\t%s\n"
                    txt += bytes(this.seq, "utf-8") + b"\n"
                    txt += b"+\n"
                    txt += bytes(this.qual,"utf-8") + b"\n"

                # Here, we must be careful as to keep the pairs. So if R1 is mapped
                # but R2 is unmapped (or the inverse), then the pair is mapped
                if this.is_read1:
                    if this.is_unmapped and this.mate_is_unmapped:
                        R1_unmapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_unmapped'] += 1
                    else:
                        R1_mapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_mapped'] += 1
                elif this.is_read2:
                    if this.is_unmapped and this.mate_is_unmapped:
                        R2_unmapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_unmapped'] += 1
                    else:
                        R2_mapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_mapped'] += 1
                else:
                    #pass
                    print("Neither R1 or R2")
                if this.is_duplicate:
                    stats['duplicated'] += 1
            pb.animate(i+1)

        print("\nNumber of entries in the BAM: %s" % str(i+1))
        R1_mapped.close()
        R2_mapped.close()
        R1_unmapped.close()
        R2_unmapped.close()
        return stats
    elif mode == "se":
        raise NotImplementedError
    else:
        raise ValueError("mode must be se or pe")


