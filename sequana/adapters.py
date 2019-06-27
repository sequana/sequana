# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Utilities to manipulate adapters

Adapters removal can be performed by many different tools such as CutAdapt,
AlienTrimmer, Trimmomatic. Unfortunately, they tend to use different formats
from FASTA to text files. Moreover outputs are generally also reported in
different formats.

Tools to extract specific adapters from FASTA files would also be handy. For
instance, you may have all your adapters in a single file.

In this module, we provide:

- tools to manipulate adapters stored in Fasta format (:class:`AdapterReader`).
- tools to export Fasta files with adapter content into other formats required
  by various adapter removal software
- A tool used to extract adapters from a FASTA file given their identifier, or
  sequence :class:`FindAdaptersFromDesign`.

Our convention is to store list of adapters in FASTA format, which can be read using
the :class:`AdapterReader`::

    from sequana import sequana_data, AdapterReader
    filename = sequana_data("adapters_Nextera_fwd.fa")
    ar = AdapterReader(filename)
    ar.get_adapter_by_index("N501")

Given a design file (see mod:`sequana.expdesign`), and a name for the type of
adapters, one can easily extract the subset of relevant adapters to be used for
a sample. Currently, the following set of adapters/design are available:

    - Nextera single and double indexing
    - Rubicon single indexing
    - PCRFree single indexing
    - TruSeq
    - NEBNext single and double indexing

Note that TruSeq index 17, 24, and 26 are missing. This is normal. Those are
"reserved" Illumina index.

For instance given a design file that gives the mapping between samples and a
set of Nextera adapters, one would use:

.. doctest::

    >>> from sequana import *
    >>> filename = sequana_data("test_expdesign_hiseq.csv")
    >>> design = ExpDesignAdapter(filename)
    >>> fa = FindAdaptersFromDesign(design, "PCRFree")
    >>> print(fa.sample_names[0])
    '553-iH2-1'
    >>> fa.get_adapters_from_sample("553-iH2-1")

See :class:`FindAdaptersFromDesign` for details.

"""
import os

from sequana import logger
from sequana.fasta import FastA
from sequana.datatools import sequana_data

logger.name = __name__

import pysam


def fasta_fwd_rev_to_columns(file1, file2=None, output_filename=None):
    """From 2 FASTA files (reverse and forward) adapters, returns 2-columns file

    This is useful for some tools related to adapter removal that takes as input
    this kind of format

    :param str filename1: FASTA format
    :param stsr filename2: FASTA format (optional)

    The files must have a one-to-one mapping
    """
    f1 = pysam.FastxFile(file1)
    if output_filename is not None:
        fout = open(output_filename, "w")

    if file2:
        f2 = pysam.FastxFile(file2)
        for read1, read2 in zip(f1, f2):
            txt = "%s %s" % (read1.sequence, read2.sequence)
            if output_filename is None:
                print(txt)
            else:
                fout.write(txt+"\n")
    else:
        for read1 in f1:
            txt = "%s" % read1.sequence
            if output_filename is None:
                print(read1.sequence)
            else:
                fout.write(txt+"\n")
    if output_filename is not None:
        fout.close()


def adapter_removal_parser(filename):
    """Parses output of AdapterRemoval software

    .. doctest::

        >>> from sequana import adapters, sequana_data
        >>> data = sequana_data("test_adapter_removal_output.txt", "testing")
        >>> results = adapters.adapter_removal_parser(data)
        >>> results["adapter1"]
        'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'

    """
    results = {}

    with open(filename, "r") as fin:
        lines = fin.readlines()
        for line in lines:
            if line.startswith("  --adapter"):
                lhs, rhs = line.split(":")
                name = lhs.strip().replace("-", "")
                sequence = rhs.strip()
                results[name] = sequence
    return results

def _get_registered_adapters():
    filenames = sequana_data('*', 'data/adapters')
    filenames = [x for x in filenames if "_fwd" in x]
    registered = [x.lstrip("adapters_").replace("_fwd.fa", "") for x in filenames]
    registered = set(registered)
    return registered


def get_sequana_adapters(type_, direction):
    """Return path to a list of adapters in FASTA format

    :param tag: PCRFree, Rubicon, Nextera
    :param type_: fwd, rev, revcomp
    :return: path to the adapter filename

    """
    # search possible types
    registered = _get_registered_adapters()
    if type_ not in registered:
        logger.error("This adapter type (%s) is not valid" % type_)
        logger.error("choose one in %s types" % registered)
        raise ValueError

    directions = ["fwd", "rev", "revcomp"]
    if direction not in directions:
        logger.error("This kind of tag (%s) is not valid" % direction)
        logger.error("choose one in %s " % directions)
        raise ValueError
    try:
        this = sequana_data("adapters_%s_%s.fa" % (type_, direction))
        logger.warning("Rename {} (remove the adapters_ prefix)".format(this))
        return this
    except:
        return sequana_data("%s_%s.fa" % (type_, direction))


class Adapter(object):
    """Class to store one adapter

    An adapter is just a sequence from a FASTA file. It contains
    an identifier, a sequence and possibly a comment.

    .. warning:: The identifier provided must not contain the starting
        ">" character, which is added automatically when needed.

    One can check if an adapter is equal to another. Only the sequence is
    checked for equality though.

    Some Sequana notation have been added in the identifier to ease retrieval
    of index's name and index's sequence::

        >NextFlex_PCR_Free_adapter1|name:1|seq:CGATGT

    Of course the string CGATGT must be found in the sequence itself.

    ::

        ar = AdapterReader(sequana_data("adapters_PCRFree_fwd.fa"))
        adapter = Adapter(ar[0])
        adapter.identifier
        adapter.comment
        adapter.index_sequence
        adapter.sequence
        adapter.name

    """
    def __init__(self, identifier, sequence=None, comment=None):
        if isinstance(identifier, dict):
            self._identifier = identifier['identifier'].strip()
            self._comment = identifier['comment'].strip()
            self._sequence = identifier['sequence'].strip()
        elif isinstance(identifier, Adapter):
            self._identifier = identifier.identifier
            self._comment = identifier.comment
            self._sequence = identifier.sequence
        else:
            if comment is None: comment = ""
            if sequence is None: sequence = ""
            self._comment = comment.strip()
            self._sequence = sequence.strip()
            self._identifier = identifier.strip()

        if self.identifier.startswith(">"):
            error = "identifier must be a string without starting > character"
            raise ValueError(error)

    def _get_sequence(self):
        return self._sequence
    def _set_sequence(self, sequence):
        self._sequence = sequence
    sequence = property(_get_sequence, _set_sequence,
        doc="R/W adapter's sequence")

    def _get_identifier(self):
        return self._identifier
    def _set_identifier(self, identifier):
        self._identifier = identifier
    identifier = property(_get_identifier, _set_identifier,
        doc="R/W adapter's identifier")

    def _get_value_from_identifier(self, tag):
        assert tag in ["name", "seq"]
        entries = [this for this in self.identifier.split('|')
                    if this.startswith(tag + ":")]
        assert len(entries) <= 1
        if len(entries) == 0:
            return None
        else:
            name = entries[0]
            return name.split(":")[1]

    def _get_name(self):
        return self._get_value_from_identifier("name")
    name = property(_get_name, doc="Read only access to the inedx name")

    def _get_seq(self):
        return self._get_value_from_identifier("seq")
    index_sequence = property(_get_seq,
        doc="Read only access to the index sequence")

    def _get_comment(self):
        return self._comment
    def _set_comment(self, comment):
        self._comment = comment
    comment = property(_get_comment, _set_comment,
        doc="R/W adapter's identifier")

    def __str__(self):
        txt = ">%(identifier)s\t%(comment)s\n"
        txt+= "%(sequence)s"

        txt = txt % {"identifier":self.identifier,
                "comment":self.comment, "sequence":self.sequence}
        return txt

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self._comment != other.comment:
            return False
        if self._identifier != other.identifier:
            return False
        if self._sequence.upper() != other.sequence.upper():
            return False
        return True


class AdapterReader(object):
    """Reader of FASTA file dedicated to adapters

    A Fasta is just a set of this kind of paired-lines::

        >Nextera_index_N501|name:N501|seq:ACGT optional comment
        ACGTACGTACGT

    where the *optional comment* is separated from the identifier by a
    tabulation.

    In the FASTA identifier, the first pipe delimits the official name (left
    hand side) from the name tag. The information on this example may
    be redundant but the *name* will be used throughout the
    **Sequana** code to ensure reproductibility.

    .. note:: sequences are all in big caps.

    .. note:: the universal adapter has no index so does not need to have the
        any tags for the name of index sequence. However, it must be
        called *Universal_Adapter*

    .. doctest::

        >>> from sequana import sequana_data, AdapterReader
        >>> filename = sequana_data("adapters_Nextera_fwd.fa")
        >>> ar = AdapterReader(filename)
        >>> candidate = ar.get_adapter_by_index_name("S505")
        >>> print(candidate[0])
        >Nextera_index_S505|name:S505|seq:GTAAGGAG
        AATGATACGGCGACCACCGAGATCTACACGTAAGGAGTCGTCGGCAGCGTC
        >>> len(ar)
        56

    .. note:: Checks for uniqueness of the identifers. It not unique, an error
        is raised

    :sources: document illumina #1000000002694  v01
    :sources: For NextFlex PCR-Free adapters, there are 48 barcodes.
        http://www.biooscientific.com/Portals/0/IEM/Bioo-Scientific-PCR-Free-Barcode-Indices-v1-1-15.pdf

    """
    def __init__(self, filename):
        """.. rubric:: Constructor

        :param str filename: the input FASTA file
        """
        if isinstance(filename, str):
            # this is not large files so we load in memory all sequences/names
            # and comments once for all. This has also the adavantage that
            # data can now be changed on the fly
            fasta = FastA(filename)
            self._data = [self._to_read(this) for this in fasta]
        elif isinstance(filename, AdapterReader):
            self._data = [self._to_read(this) for this in filename._data]
        elif isinstance(filename, list):
            self._data = [self._to_read(this) for this in filename]

        self._sanity_check()

    def __len__(self):
        return len(self._data)

    def __getitem__(self, i):
        return self._data[i]

    def __repr__(self):
        txt = 'AdapterReader. List of %s adapters' % len(self)
        return txt

    def _get_field(self, tag):
        # usage: _get_field("seq")
        matches = []
        tag = tag + ":"
        for this in self.identifiers:
            if tag in this:
                for field in this.split("|"):
                    field = field.strip()
                    if field.startswith(tag):
                        matches.append(field.replace(tag, ""))
            else:
                matches.append(None)
        return matches

    def _get_identifiers(self):
        return [this.identifier for this in self._data]
    identifiers = property(_get_identifiers)

    def _get_seq(self):
        return [this.sequence for this in self._data]
    sequences = property(_get_seq)

    def _get_comments(self):
        return [this.comment for this in self._data]
    comments = property(_get_comments)

    def _get_index_sequences(self):
        return self._get_field("seq")
    index_sequences = property(_get_index_sequences)

    def _get_index_names(self):
        return self._get_field("name")
    index_names = property(_get_index_names)

    def _sanity_check(self):
        """Check that all identifiers are unique"""
        if len(set(self.identifiers)) != len(self.identifiers):
            import collections
            identifiers = [k for k,v in collections.Counter(self.identifiers).items() if v>1]
            msg = "Found identical identifiers in fasta sequences\n"
            msg += "Check those identifiers for duplicates %s " % identifiers
            raise ValueError(msg)

    def get_adapter_by_sequence(self, subsequence):
        """Return one or several adapters with sub-sequence in their sequence

        :param str subsequence: a string (ACGT letters)
        :return: name and sequence in FASTA format that have the user *sequence*
            contained in their sequence

        If the subsequence is short, it may return more than 1 adapters. Besides,
        the sequence is searched for without position information right now.
        """
        found = [i for i,x in enumerate(self.sequences) if subsequence in x]
        if len(found) == 0:
            return None

        adapters = []
        for index in found:
            this = self._data[index]
            this_adapter = Adapter(identifier=this.identifier,
                                   sequence=this.sequence,
                                   comment=this.comment)
            adapters.append(this_adapter)
        return adapters

    def get_adapter_by_identifier(self, text):
        """Return adapter whose identifier matches the user text

        :param index_identifier: the unique index identifier to be found.
            If several sequence do match, this is an error meaning the fasta
            file with all adapters is not correctly formatted.
        :return: the adapter that match the index_name (if any) otherwise
            returns None

        """
        adapters = []
        for this in self._data:
            if text == this.identifier or text in this.identifier.split("|"):
                this_adapter = Adapter(identifier=this.identifier,
                                       sequence=this.sequence,
                                       comment=this.comment)
                adapters.append(this_adapter)

        if len(adapters) == 0:
            raise ValueError("No matching identifier found with this pattern: %s" % text)
        elif len(adapters) == 1:
            return adapters[0]
        else:
            raise ValueError("Found two adapters matching the identifier. This should never happen")

    def _get_adapter_by_index(self, index_name, prefix):
        """Return adapter corresponding to the unique index

        :param index_name: the unique index name to be found. If several
            sequence do match, this is an error meaning the fasta file
            with all adapters is not correctly formatted.
        :return: an instance of :class:`Adapter` if index_name match an
            adapter; returns None otherwise

        ::

            from sequana import sequana_data, AdapterReader
            filename = sequana_data("adapters_Nextera_fwd.fa")
            ar = AdapterReader(filename)
            ar.get_adapter_by_identifier("N712")

        """
        # there should be only one
        adapters = []
        for this in self._data:
            if prefix + str(index_name) in this.identifier.split("|"):
                this_adapter = Adapter(identifier=this.identifier,
                                       sequence=this.sequence,
                                       comment=this.comment)
                adapters.append(this_adapter)

        if len(adapters) == 0:
            return None
        elif len(adapters) >= 2:
            logger.warning("Found two adapters matching index {}. This may happen with Nextera adapters".format(index_name))
        return adapters

    def get_adapter_by_index_name(self, index_name):
        """Return adapter for the index name provided

        Can be used only if the identifier contains the tag::

            |name:an_index_to_be_found

        For instance::

            >Nextera_blabal|name:N505|seq:ACGT
            >Nextera_blabal|seq:ACGT|name:N505
            >Nextera_blabal|name:N505

        are valid identifiers

        """
        return self._get_adapter_by_index(index_name, prefix="name:")

    def get_adapter_by_index_seq(self, index_name):
        """See :meth:`get_adapter_by_index_name`."""
        return self._get_adapter_by_index(index_name, prefix="seq:")

    def _to_read(self, this):
        from easydev import AttrDict
        d = AttrDict()
        d.sequence = this.sequence
        if this.comment is None:
            this.comment  = ""

        d.comment = this.comment
        try:
            #pysam format
            d.identifier = this.name
        except:
            # this class convention
            d.identifier = this.identifier
        return d

    def to_dict(self):
        """Returns dictionary with key as identifier and values as
        list with comments and sequences
        """
        d1 = [(this.identifier, [this.comment, this.sequence])
              for this in self._data]
        return dict(d1)

    def __eq__(self, other):
        other = AdapterReader(other)
        d1 = self.to_dict()
        d2 = other.to_dict()
        return d1 == d2

    """def __reverse_comp(self, this):
        from sequana.sequence import DNA
        if this.startswith("seq:"):
            tag, seq = this.split(":")
            return tag + ":" + DNA(seq[:]).get_reverse_complement()
        else:
            return this

    def __reverse(self, this):
        if this.startswith("seq:"):
            tag, seq = this.split(":")
            return tag + ":" + seq[::-1]
        else:
            return this
    """
    def reverse(self):
        """Reverse all sequences inplace

        .. doctest::

            >>> from sequana import sequana_data, AdapterReader
            >>> filename = sequana_data("adapters_Nextera_fwd.fa")
            >>> filename2 = sequana_data("adapters_Nextera_rev.fa")
            >>> ar = AdapterReader(filename)
            >>> ar2 = AdapterReader(filename2)
            >>> ar.reverse()
            >>> ar == ar2
            True

        """
        # Reverse in place
        for this in self._data:
            this.sequence = this.sequence[::-1]
            fields = this.identifier.split("|")
            # change API version 0.7.2
            # we now store the original sequence instead of reverse
            # identifier = "|".join([self._reverse(field) for field in fields])
            identifier = "|".join([field for field in fields])
            this.identifier = identifier

    def reverse_complement(self):
        """Reverse-complement all sequences inplace

        ::

            >>> from sequana import sequana_data, AdapterReader
            >>> filename = sequana_data("adapters_Nextera_fwd.fa")
            >>> filename = sequana_data("adapters_Nextera_revcomp.fa")
            >>> ar = AdapterReader(filename)
            >>> ar.reverse_complement()
            >>> ar.to_fasta()
            >>> ar == ar2

        """
        # no need for a fast implementation (less than 1000 short reads)
        from sequana.sequence import DNA
        for this in self._data:
            this.sequence = DNA(this.sequence[:]).get_reverse_complement()
            fields = this.identifier.split("|")
            # change API version 0.7.2
            # we now store the original sequence instead of reverse
            # identifier = "|".join([self._reverse_comp(field) for field in fields])
            identifier = "|".join([field for field in fields])
            this.identifier = identifier

    def to_fasta(self, filename):
        """Save sequences into fasta file"""
        with open(filename, "w") as fout:
            for i, this in enumerate(self._data):
                if i > 0:
                    fout.write("\n")
                comment = this.comment
                if comment is None or comment == "":
                    fout.write(">%s\n%s" % (this.identifier, this.sequence))
                else:
                    fout.write(">%s\t%s\n%s" % (this.identifier, comment,
                                            this.sequence))


class FindAdaptersFromDesign(object):
    """Extract adapter(s) corresponding to an experimental design file

    Used by sequana main script to build the adapter files for multi-samples
    projects as input to various adapter removal software.

    """
    def __init__(self, design_filename, adapters):
        """.. rubric:: Constructor

        :param str design_filename: a CSV file that is compatible
            with our :class:`sequana.expdesign.ExpDesignAdapter`
        :param adapters: the type of adapters (PCRFree, Nextera,
            Rubicon, TruSeq, SMARTer, Small)

        The files of adapters are stored in Sequana and accessible with the
        sequana_data function. So, for instance if adapters is set to Nextera,
        the following file is used to identify the adapters::

            sequana_data("adapters_Nextera_fwd.fa")

        New adapters files can be added on request. See resources/data/adapters
        for the full list. You can also use::

            from sequana.adapters import _get_registered_adapters
            _get_registered_adapters()
        """
        from sequana.expdesign import ExpDesignAdapter
        self.design = ExpDesignAdapter(design_filename)

        if self.design.df.index.name == "Sample_ID" or \
            "Sample_ID" in self.design.df.columns:
            self.design.df.set_index("Sample_ID", inplace=True)
        else:
            raise ValueError("Incorrect design file. Missing Sample_ID field")

        self.adapters = adapters

        try:
            file1 = sequana_data("adapters_%s_fwd.fa" % adapters)
            logger.warning("rename your file removing prefix adatper")
        except:
            file1 = sequana_data("%s_fwd.fa" % adapters)

        try:
            file2 = sequana_data("adapters_%s_revcomp.fa" % adapters)
            logger.warning("rename your file removing prefix adatper")
        except:
            file2 = sequana_data("%s_revcomp.fa" % adapters)

        self._adapters_fwd = AdapterReader(file1)
        self._adapters_revc = AdapterReader(file2)  # !!! revcomp

    def _get_samples(self):
        return list(self.design.df.index)
    sample_names = property(_get_samples,
        doc="return all sample names contained in the design file")

    def get_sample(self, sample_name):
        """Return basic info about the sample name (from the design file)"""
        if sample_name not in self.design.df.index:
            raise ValueError("%s not valid. Use one of %s" % (sample_name,
                                                              self.sample_names))

        data = self.design.df.loc[sample_name]
        if data.ndim == 1: # the expected pandas.Series
            return data
        else:
            # Check that we have duplicates
            # Indeed, for HiSeq design with Index1_Seq and Index2_Seq, it may happen
            # that the sample names are duplicated (RNA-seq exp) on  different Lanes.
            checkme = data.drop_duplicates(["SampleRef", "Index1_Seq", "Index2_Seq"])
            if len(checkme) > 1:
                raise ValueError("Found several instance of sample " +
                        "name {} in the design file".format(sample_name))
            # return only the first instance
            return data.iloc[0]

    def get_adapters_from_sample(self, sample_name):
        """Return a dictionary with adapters corresponding to the sample name

        :param str sample_name: a valid sample name as found in the design
            file. One can check the content of the :attr:`sample_names`
            attribute.
        :return: a dictionary with the adapters in forward and reverse
            complement for index1 and index2 (if relevant).

        ::

            >Universal_Adapter|name:universal
            AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
            >Nextera_transposase_seq_1|name:transposase_seq_1
            TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
            >Nextera_transposase_seq_2|name:transposase_seq_2
            GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
            >Nextera_index_N501|name:N501|seq:TAGATCGC

        If sample name is given, it will figure out that for instance the index
        N501 is required. In addition, all sequences with identifier without the
        |seq: tag in their name will be added; So, here we will also have the
        universal and the two transposases


        """
        data = self.get_sample(sample_name)

        res = {'index1': {}, 'index2': {}}

        # Index1_Seq must always be present. This is part of the API of the
        # ExpDesignAdapter class. However, Index2_ID may not always be present
        # In which case index2 remains empty

        # Then, two types of design are accepted, using the adapter index
        # ID or the sequence itself. The sequence is more robust since
        # experimentalist may change the ID (but not the seq). So we start with
        # the sequence first.

        for column in ["Index1_Seq", "Index2_Seq"]:
            if column not in data.index:
                continue
            key = column.split("_")[0].lower()
            index = data.loc[column]
            if index is None:
                continue
            seq = self._adapters_fwd.get_adapter_by_index_seq(index)
            if seq is None:
                raise ValueError("Found no index for %s (sample %s)" % (index,sample_name))

            # drop potential duplicates
            uniques = set([this.sequence for this in seq])
            if len(uniques) > 1:
                raise ValueError("Found 2 sequences with index %s " % (index, sample_name))

            # If there are duplicates, or not, just take first element
            res[key]['fwd'] = seq[0]

            # Reverse comp version
            seq = self._adapters_revc.get_adapter_by_index_seq(index)
            if seq is None: #Index2 is optional so no error raised
                pass
            else:
                uniques = set([this.sequence for this in seq])
                if len(uniques) > 1:
                    raise ValueError("Found 2 sequences with index %s " % (index, sample_name))
                res[key]['revc'] = seq[0]

        # If Index1_Seq not in the index, then we should use the IDs
        if "Index1_Seq" not in data.index:
            if "Index1_ID" in data.index:
                logger.warning("Usage of IDs for indexing adapters is deprecated. See adapters module")
                index1 = data.loc['Index1_ID']
                res['index1']['fwd'] = self._adapters_fwd.get_adapter_by_index_name(index1)
                res['index1']['revc'] = self._adapters_revc.get_adapter_by_index_name(index1)
            if "Index2_ID" in data.index:
                logger.warning("Usage of IDs for indexing adapters is deprecated. See adapters module")
                index2 = data.loc['Index2_ID']
                res['index2']['fwd'] = self._adapters_fwd.get_adapter_by_index_name(index2)
                res['index2']['revc'] = self._adapters_revc.get_adapter_by_index_name(index2)

        # adapters to be found but not a barcode, whih have the "|seq:" string
        # in it
        for name, identifier in zip(self._adapters_fwd.index_names,
                                    self._adapters_fwd.identifiers):
            if "|seq:" in identifier:
                pass
            else:
                res[name] = {}
                res[name]['fwd'] = \
                    self._adapters_fwd.get_adapter_by_index_name(name)[0]
                res[name]['revc'] = \
                    self._adapters_revc.get_adapter_by_index_name(name)[0]

        # FIXME changes the dictionary in the loop. May not be wise
        res = dict([(k,v) for k,v in res.items() if len(v)!=0])
        return res

    def check(self):
        found = 0
        for sample in self.sample_names:
            try:
                self.get_adapters_from_sample(sample)
                found += 1
            except:
                logger.error("No index found for sample %s" % sample)
        if found == 0:
            raise ValueError("None of the sample match any of the adapters")

    def save_adapters_to_fasta(self, sample_name, output_dir='.'):
        """Get index1, index2 and other standard adapters"""
        adapters = self.get_adapters_from_sample(sample_name)

        file_fwd = output_dir + os.sep + "%s_adapters_fwd.fa"% sample_name
        with open(file_fwd, "w") as fout:
            keys = sorted(adapters.keys())
            for k in keys:
                fout.write(str(adapters[k]['fwd'])+"\n")

        file_revc = output_dir + os.sep + "%s_adapters_revcomp.fa" % sample_name
        with open(file_revc, "w") as fout:
            keys = sorted(adapters.keys())
            for k in keys:
                fout.write(str(adapters[k]['revc'])+"\n")


        return file_fwd, file_revc
