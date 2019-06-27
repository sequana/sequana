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
"""Module to handle experimental design files (adapters)


Sequencers or experimentalists create so-called design files to
store information about the sequencing experiment. For example the
name of the samples, the sample well, and the index sequences.

The format used to store the design information may vary from one sequencing
platform to another. The design file can be used for many different purposes.
Currently, we only use them to retrieve information about adapters.

Since there are different formats possible, we decided on a minimal and common
set of information. For now, this is a CSV file with the following minimal header::

    Sample_ID, Index1_Seq

or, for double-indexing::

    Sample_ID, Index1_Seq, Index2_Seq


Users should only use the :class:`ExpDesignAdapter` class, which understands
the different formats. Currently, the design file created by MiSeq Illumina
machine


"""
import shlex
import io

from sequana.lazy import pandas as pd
from sequana import logger

logger.name = __name__


__all__ = ["ExpDesignAdapter", "ExpDesignMiSeq", "ExpDesignHiSeq",
           "ExpDesignBase"]


class ExpDesignAdapter(object):
    """Generic Experimental design class

    This class is used to store the mapping between sample ID and adapter
    used.

    The design information is stored as a dataframe in the attribute
    :attr:`df`.

    The input format is not unique. There are currently 3 different inputs
    possible as defined in

    - :class:`ExpDesignGeneric`
    - :class:`ExpDesignMiSeq`
    - :class:`ExpDesignHiSeq` (2500)

    The dataframe index is the list of sample identifiers (Sample_ID).
    The columns contain at least the following::

        Index1_Seq, Index1_ID, Index2_Seq, Index2_ID

    Example::

        from sequana import *
        filename = sequana_data('test_test_expdesign_hiseq.csv')
        eda = ExpDesignAdapter(filename)

    """
    def __init__(self, filename, verbose=True):
        """.. rubric:: constructor

        :param str filename: the input design file. Can also be an instance of
            :class:`ExpDesignAdapter` itself.
        :param bool verbose:

        """
        self.verbose = verbose
        try:
            self.df = filename.df.copy()
            self.name = "unset"
        except:
            self._factory(filename)

    def _factory(self, filename):
        try:
            exp = ExpDesignMiSeq(filename)
            if self.verbose:
                logger.info("Found MiSeq design file")
        except Exception as err:
            try:
                exp = ExpDesignHiSeq(filename)
                if self.verbose:
                    logger.info("Found HiSeq design file")
            except:
                try:
                    exp = ExpDesignGeneric(filename)
                    if self.verbose:
                        logger.info("Found Generic Sequencer design file")
                except:
                    msg = "Input file could not be read or interpreted"
                    raise IOError(msg)
        self.df = exp.df
        self.name = exp.name

    def __repr__(self):
        txt = "ExpDesignAdapter: %s experiments\n" % len(self.df)
        txt += " read with %s\n" % self.name
        try: txt += " Number of lanes: %s \n" % len(self.df.Lane.unique())
        except: pass
        try: txt += " Number of sample reference: %s" % len(self.df.SampleRef.unique())
        except: pass
        return txt


class ExpDesignBase(object):
    """The Base class for all ExpDesignAdapter classes

    The input filename must be a CSV file with at least the following column in
    the header::

        Sample_ID

    Derived class must define at least **Index1_Seq**  and possibly **Index2_Seq**.

    Examples of specialised classes are :class:`ExpDesignMiSeq`,
    :class:`ExpDesignHiSeq`.

    """
    def __init__(self, filename, sep=","):
        self.sep = sep
        self.filename = filename
        self.adapter_type = "unset"
        self.name = "ExpDesignBase"
        self.instrument = "undefined"

    def check(self):
        """Check the presence of the Sample_ID column"""
        for this in ["Sample_ID"]:
            if this not in self.df.columns:
                print(self.df.columns)
                raise KeyError("%s not found. " % this)

    def read(self):
        """Read a CSV file"""
        self.df = pd.read_csv(self.filename, sep=self.sep)

    def __repr__(self):
        txt = "%s: %s entries\n" % (self.name, len(self.df))
        txt += "adapter type: %s\n" % (self.adapter_type)
        txt += "Instrument Type: %s " % (self.instrument)
        return txt


class ExpDesignGeneric(ExpDesignBase):
    """Generic experimental design format

    ::

        Sample_ID, Index1_Seq, Index2_Seq

    :: 

    """
    def __init__(self, filename):
        super(ExpDesignGeneric, self).__init__(filename)
        self.name = "ExpDesignGeneric"

        self.read()
        self.df.rename(columns={"index1": "Index1_ID", "index2": "Index2_ID",
            "sample_name": "Sample_ID"},  inplace=True)
        self.check()


class ExpDesignHiSeq(ExpDesignBase):
    """Dedicated experimental design format created by a demultiplexing soft.

    This format is used by a demultiplex software used locally at biomics
    platform. The format of the header is::

        FCID,Lane,SampleID,SampleRef,Index Seq, Description,Control,Recipe,Operator,Project

    This is a format that may change in the future.

    The SampleID is converted into Sample_ID, "Index Seq". Note that "Index Seq"
    may be empty, or filled with an index sequence, or 2 index sequences
    separated by a "-" sign.

    note also FCID = flowcell ID

    """
    def __init__(self, filename, sep=","):
        super(ExpDesignHiSeq, self).__init__(filename, sep=sep)
        self.name = "ExpDesignHiSeq"
        self.read()

        for this in ["FCID","Lane","SampleID","Index Seq"]:
            if this not in self.df.columns:
                raise ValueError("Invalid header. %s column is missing" % this)

        self.df.rename(columns={"SampleID":"Sample_ID", "Index Seq":"Index_Seq"},
                       inplace=True)

        index1 = []
        index2 = []
        for this in self.df.Index_Seq.values:
            if isinstance(this, str):
                indices = this.split("-")
                index1.append(indices[0])
                if len(indices) == 2:
                    index2.append(indices[1])
                else:
                    index2.append(None)
            else:
                index1.append(None)
                index2.append(None)
        self.df['Index1_Seq'] = index1
        self.df['Index2_Seq'] = index2
        self.df.drop("Index_Seq", axis=1, inplace=True)

        self.check()


class ExpDesignIllumina(ExpDesignBase):
    def __init__(self, filename):
        super(ExpDesignIllumina, self).__init__(filename)

    def scanner(self):
        data = {}
        # shlex removes all white lines and split by return carriage
        # strip is also applied
        rawdata = shlex.split(open(self.filename, "r"))
        for line in rawdata:
            # sometimes, IEM will store the ;;; at the end
            # so we can get [HEADER];;;;;;;;;;;
            if line.startswith('[') and "]" in line:
                line = line.strip(";").strip(",").strip()
                currentkey = line.replace("[", "").replace("]", "")
                data[currentkey] = []
            else:
                data[currentkey].append(line)

        for key in data.keys():
            data[key] = "\n".join(data[key])

        for this in ["Header", "Reads", "Settings", "Data"]:
            if this not in data.keys():
                logger.warning("%s not found in the DesignExpMiSeq file" % this)

        self.data = data
        self.df = pd.read_csv(io.StringIO(data["Data"]))

        ncols = [8, 9,10, 12]
        if self.df.shape[1] not in ncols:
            self.df = pd.read_csv(io.StringIO(data["Data"]), ";")
            if self.df.shape[1] not in ncols:
                logger.warning("Data section must have 10 or 12 columns. Check the samplesheet")

        # Fixes https://github.com/sequana/sequana/issues/507
        self.df["Sample_ID"] = self.df["Sample_ID"].astype(str)

        self.df.rename(columns={"I7_Index_ID":"Index1_ID", "index":"Index1_Seq",
            "I5_Index_ID": "Index2_ID", "index2":"Index2_Seq"},
                       inplace=True)

    def to_demultiplex(self, filename, FCID, double_indexing=True):

        # FCID,Lane,SampleID,SampleRef,Index Seq,Description,Control,Recipe,Operator,Project
        # HYVN3BCXY,1,177-SQDC,170629_C-dHumieres,CTCTCTAC-GCGTAAGA,,N,V2.0,PF1,170629_C-dHumieres

        if "Lane" not in self.df.columns:
            self.df['Lane'] = 1

        if double_indexing is False:
            raise NotImplementedError

        with open(filename, "w") as fout:
            fout.write("FCID,Lane,SampleID,SampleRef,Index Seq,Description,Control,Recipe,Operator,Project\n")
            for this in self.df.iterrows():
                data = this[1]
                txt = "{},{},{},{},{}-{},,N,,PF1,{}\n".format(
                        FCID,
                        data.Lane,
                        data.Sample_ID,
                        data.Sample_Project,
                        data.Index1_Seq,
                        data.Index2_Seq, 
                        data.Sample_Project)
                fout.write(txt)



class ExpDesignHiSeq2500(ExpDesignIllumina):
    def __init__(self, filename):
        super(ExpDesignHiSeq2500, self).__init__(filename)
        self.name = "ExpDesignHiSeq2500"

        self.scanner()


        # The name of the Index_ID is not standard....
        # Depends on the experimentalist because a prefix may be added.
        # One known prefix is NF. We agreed that future prefix must end with an
        # underscore so that it can be removed. Since ID may contain letters
        # (e.g.S501), it would be impossible otherwise to split the prefix from
        # the index.
        self.df["Index1_ID"] = self.df["Index1_ID"].apply(
                lambda x: x.replace("NF", ""))
        self.df["Index1_ID"] = self.df["Index1_ID"].apply(
                lambda x: x.split("_",1)[-1])
        try:
            self.df["Index1_ID"] = self.df["Index1_ID"].astype(int)
        except:
            pass

        if "Index2_ID" in self.df.columns:
            self.df["Index2_ID"] = self.df["Index2_ID"].apply(
                lambda x: x.replace("NF", ""))
            self.df["Index2_ID"] = self.df["Index2_ID"].apply(
                lambda x: x.split("_",1)[-1])
            try:
                self.df["Index2_ID"] = self.df["Index2_ID"].astype(int)
            except:
                pass

        # Figure out the type of adapters if possible
        try:
            header = self.data['Header']
            assay = [x for x in header.split('\n') if x.startswith("Assay")]
            assay = assay[0]
            items = assay.split(',')
            if items[0] == "Assay":
                self.adapter_type = assay.split(",")[1]
            elif items[0] == "Instrument Type":
                self.instrument = assay.split(",")[1]
            else:
                items = assay.split(';')
                self.adapter_type = assay.split(";")[1]
        except:
            pass

        self.check()


class ExpDesignMiSeq(ExpDesignIllumina):
    """Dedicated experimental design format from Illumina MiSeq sequencers

    This MiSeq design format has the following format::

        [Header]
        blabla

        [Reads]
        blabla

        [Settings]
        blabla

        [Data]
        blabla


    In the Data section, a CSV file is to be found with the following header::

        Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,
        Sample_Project,Description

    The index column may be prefixed. For instance as
    NFXX where XX is the index so NF should be dropped.

    If double-indexing, the header is::

        Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,
        index,I5_Index_ID,index2,Sample_Project,Description

    ::

        filename = sequana_data("test_expdesign_miseq_illumina_1.csv")
        ff = ExpDesignMiSeq(filename)
        ff.df


    """
    def __init__(self, filename):
        super(ExpDesignIllumina, self).__init__(filename)
        self.name = "ExpDesignMiSeq"

        self.scanner()



        # The name of the Index_ID is not standard....
        # Depends on the experimentalist because a prefix may be added.
        # One known prefix is NF. We agreed that future prefix must end with an
        # underscore so that it can be removed. Since ID may contain letters
        # (e.g.S501), it would be impossible otherwise to split the prefix from
        # the index.
        self.df["Index1_ID"] = self.df["Index1_ID"].apply(
                lambda x: x.replace("NF", ""))
        self.df["Index1_ID"] = self.df["Index1_ID"].apply(
                lambda x: x.split("_",1)[-1])
        try:
            self.df["Index1_ID"] = self.df["Index1_ID"].astype(int)
        except:
            pass

        if "Index2_ID" in self.df.columns:
            self.df["Index2_ID"] = self.df["Index2_ID"].apply(
                lambda x: x.replace("NF", ""))
            self.df["Index2_ID"] = self.df["Index2_ID"].apply(
                lambda x: x.split("_",1)[-1])
            try:
                self.df["Index2_ID"] = self.df["Index2_ID"].astype(int)
            except:
                pass

        # Figure out the type of adapters if possible
        try:
            header = self.data['Header']
            assay = [x for x in header.split('\n') if x.startswith("Assay")]
            assay = assay[0]
            items = assay.split(',')
            if items[0] == "Assay":
                self.adapter_type = assay.split(",")[1]
            else:
                items = assay.split(';')
                self.adapter_type = assay.split(";")[1]
        except:
            pass

        self.check()
