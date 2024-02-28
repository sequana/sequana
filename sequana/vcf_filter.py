#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
# CODE KEPT FOR NOW FOR BOOK KEEPING
# SHOULD BE MERGED INTO freebayes_vcf_filter


class VCF_mpileup_4dot1(VCFBase):
    """VCF filter dedicated to version 4.1 and mpileup

    This filter can be used to filter VCF created with mpileup.

    By default, no filters are applied.

    One need to define the filters.

    There are 3 hard-coded filters that can be applied by settings these
    attributes to True:

    * apply_af1_filter
    * apply_dp4_filter: includes minimum_depth, minimu_depth_strand and ratio
      filters (see later)
    * apply_indel_filter


    Terminology:

    * depth        : number of reads matching the position (4)
    * depth_strand : number of reads matching the position per strand (2)
    * ratio        : ratio of first to second base call (0.75)
    * quality      : variant quality (e.g. 50)
    * map_quality  : mapping quality (e.g. 30)
    * af1          : allele frequency (you would expect an AF of 1
        for haploid SNPs). Here set to 0.95 by default
    * strand_bias  : p-value for strand bias.
    * map_bias     : p-value for mapping bias.
    * tail_bias    : p-value for tail distance bias.


    :resources: https://github.com/sanger-pathogens/vr-codebase/ (see
        Variant/Evaluator/Pseudosequence.pm file)

    """

    def __init__(self, filename, **kwargs):
        """
        Filter vcf file with a dictionnary.
        It takes a vcf file as entry.

        You can filter a tag within the INFO list using one of those syntax
        (using the DP tag as an example):

            DP<30
            DP<=30
            DP>30
            DP>=30

        For some tags, you want to keep values within or outside a range of values.
        Tou can then use the & and | characters::

            DP<30|>60  # to keep only values in the ranges [0-30] and [60-infinite]
            DP>30&<60  # to keep only values in the range [30-60]

        Some tags stores a list of values. For instance DP4 contains 4 values.
        To filter the value at position 1, use e.g.::

            DP4[2]<4

        you can use the same convention for the range as above::

            DP4[2]>4&<1000

        you may also need something like:

            sum(DP4[2]+DP4[3]) <2

        Usage example::

            from sequana import vcf_filter
            v = vcf_filter.VCF(filename)
            v.vcf.apply_indel_filter = True
            v.vcf.filter_dict['INFO']['PV4[0]'] = "<0.001"
            v.vcf.filter_dict['INFO']['PV4[2]'] = "<0.001"
            v.vcf.filter_dict['INFO']['PV4[3]'] = "<0.001"
            v.vcf.filter_dict['INFO']['MQ'] = "<30"
            v.vcf.apply_dp4_filter = True
            v.vcf.apply_af1_filter = True
            v.vcf.filter_dict['QUAL'] = 50

            v.vcf.filter_vcf("test.vcf")


        """
        super().__init__(filename, **kwargs)
        self.filter_dict = {"QUAL": -1, "INFO": {}}

        self.apply_dp4_filter = False
        self.apply_af1_filter = False
        self.apply_indel_filter = False
        self.dp4_minimum_depth = 4
        self.dp4_minimum_depth_strand = 2
        self.dp4_minimum_ratio = 0.75
        self.minimum_af1 = 0.95

    def _filter_info_field(self, info_value, threshold):
        # Filter the line if assertion info_value compare to threshold
        # is True. for instance,
        # info_value = 40 and thrshold="<30"
        # 40 is not <30 so should return False

        if "&" in threshold:
            exp1, exp2 = threshold.split("&")
            exp1 = exp1.strip()
            exp2 = exp2.strip()
            return self._filter_info_field(info_value, exp1) and self._filter_info_field(info_value, exp2)

        if "|" in threshold:
            exp1, exp2 = threshold.split("|")
            exp1 = exp1.strip()
            exp2 = exp2.strip()
            return self._filter_info_field(info_value, exp1) or self._filter_info_field(info_value, exp2)

        if threshold.startswith("<"):
            if threshold.startswith("<="):
                if info_value <= float(threshold[2:]):
                    return True
            elif info_value < float(threshold[1:]):
                return True
        elif threshold.startswith(">"):
            if threshold.startswith(">="):
                if info_value >= float(threshold[2:]):
                    return True
            elif info_value > float(threshold[1:]):
                return True
        return False

    def _get_variant_tag(self, variant):
        return "{}/{}".format(variant.CHROM, variant.POS)

    def _filter_line(self, vcf_line, filter_dict=None, iline=None):
        """

        return False if the variant should be filter

        """
        VT = self._get_variant_tag(vcf_line)

        if filter_dict is None:  # prama: no cover
            # a copy to avoid side effects
            filter_dict = self.filter_dict.copy()

        if filter_dict["QUAL"] != -1 and vcf_line.QUAL < filter_dict["QUAL"]:
            logger.debug("{} filtered variant with QUAL below {}".format(VT, filter_dict["QUAL"]))
            return False

        if self.apply_indel_filter:
            if self.is_indel(vcf_line) is True:  # pragma: no cover
                logger.debug("{}: filter out line {} (INDEL)".format(VT, iline))
                return False

        # DP4
        if self.apply_dp4_filter and "DP4" in vcf_line.INFO:
            status = self.is_valid_dp4(
                vcf_line,
                self.dp4_minimum_depth,
                self.dp4_minimum_depth_strand,
                self.dp4_minimum_ratio,
            )
            if not status:
                logger.debug("{}: filter out DP4 line {} {}".format(VT, iline, vcf_line.INFO["DP4"]))
                return False

        # AF1
        if self.apply_af1_filter and "AF1" in vcf_line.INFO:
            status = self.is_valid_af1(vcf_line, self.minimum_af1)
            if not status:
                logger.debug("{}: filter out AF1 {} on line {}".format(VT, vcf_line.INFO["AF1"], iline))
                return False

        for key, value in filter_dict["INFO"].items():
            # valid expr is e.g. sum(DP4[2],DP4[0])
            # here, we first extract the variable, then add missing [ ]
            # brackets to make a list and use eval function after setting
            # the local variable DP4 in the locals namespace
            # PV4 skip non morphic cases (no need to filter)
            if key == "PV4" and not self.is_polymorphic(vcf_line):  # pragma no cover
                return True

            # Filter such as " sum(DP[0], DP4[2])<60 "
            if key.startswith("sum("):
                # add the missing [] to create an array

                expr = key.replace("sum(", "sum([")[0:-1] + "])"

                # identify the key
                mykey = expr[5:].split("[")[0]
                lcl = locals()
                lcl[mykey] = vcf_line.INFO[mykey]
                result = eval(expr)
                if self._filter_info_field(result, value):
                    logger.debug("{} filtered variant {},{} with value {}".format(VT, result, expr, value))
                    return False
                else:
                    return True

            # key could be with an index e.g. "DP4[0]<4"
            if "[" in key:
                if "]" not in key:  # pragma: no cover
                    raise ValueError("Found invalid filter %s" % key)
                else:
                    key, index = key.split("[", 1)
                    key = key.strip()
                    index = int(index.replace("]", "").strip())
            else:
                index = 0

            # otherwise, this is probably a valid simple filter such as "DP<4"
            try:
                if type(vcf_line.INFO[key]) != list:
                    if self._filter_info_field(vcf_line.INFO[key], value):
                        val = vcf_line.INFO[key]
                        logger.debug("{}: filtered variant {},{} with value {}".format(VT, key, value, val))
                        return False
                else:
                    Nlist = len(vcf_line.INFO[key])
                    if index > Nlist - 1:  # pragma: no cover
                        raise ValueError("Index must be less than %s (starts at zero)" % Nlist)
                    if self._filter_info_field(vcf_line.INFO[key][index], value):
                        return False
            except KeyError:  # pragma no cover
                logger.debug("The information key {} doesn't exist in VCF file (line {}).".format(key, iline + 1))
        return True

    def is_polymorphic(self, variant):
        if not len(variant.ALT):
            return False

        return not str(variant.ALT[0].value).strip() == "."

    # overwrite behaviour of Variant.is_indel
    def is_indel(self, variant):
        if not len(variant.ALT):
            return False

        if "," in str(variant.ALT[0].value):  # pragma: no cover
            return True
        if "," in str(variant.REF[0]):  # pragma: no cover
            return True
        if "INDEL" in variant.INFO.keys() and variant.INFO["INDEL"]:
            return True
        return False

    def is_valid_af1(self, variant, minimum_af1=0.95):
        if "AF1" not in variant.INFO.keys():  # pragma: no cover
            return True

        if self.is_polymorphic(variant):
            if variant.INFO["AF1"] < minimum_af1:
                return False
            else:
                return True
        else:
            if variant.INFO["AF1"] > 1 - minimum_af1:
                return False
            else:
                return True

    def is_valid_dp4(self, variant, minimum_depth, minimum_depth_strand, minimum_ratio):
        """return true if valid"""

        if "DP4" not in variant.INFO.keys():  # pragma: no cover
            return True

        ref_forward = variant.INFO["DP4"][0]
        ref_reverse = variant.INFO["DP4"][1]
        alt_forward = variant.INFO["DP4"][2]
        alt_reverse = variant.INFO["DP4"][3]

        count_reference = ref_forward + ref_reverse
        count_alt = alt_forward + alt_reverse
        # set forward ratio
        forward = float(ref_forward + alt_forward)

        if forward > 0:
            ratio_forward_reference = ref_forward / forward
            ratio_forward_alt = alt_forward / forward
        else:
            ratio_forward_reference = 0
            ratio_forward_alt = 0

        # set reverse ratio
        reverse = float(ref_reverse + alt_reverse)

        if reverse > 0:
            ratio_reverse_reference = ref_reverse / reverse
            ratio_reverse_alt = alt_reverse / reverse
        else:
            ratio_reverse_reference = 0
            ratio_reverse_alt = 0
        if not self.is_polymorphic(variant):
            # dealing with non polymorphic site (i.e; VCF's ALT
            # field equals "."

            # reference depth test for the reference allele
            if count_reference < minimum_depth:  # pragma: no cover
                return False

            # forward strand depth test for the reference allele
            if ref_forward < minimum_depth_strand:
                return False

            # reverse strand depth test for the reference allele
            if ref_reverse < minimum_depth_strand:
                return False

            # forward ratio test for the reference allele
            if ratio_forward_reference < minimum_ratio:
                return False

            # reverse ratio test for the reference allele
            if ratio_reverse_reference < minimum_ratio:
                return False
        else:
            # reference depth test for the alternate allele
            if count_alt < minimum_depth:
                return False

            # forward strand depth test for the alternate allele
            if alt_forward < minimum_depth_strand:
                return False

            # reverse strand depth test for the alternate allele
            if alt_reverse < minimum_depth_strand:
                return False

            # forward ratio test for the reference allele
            if ratio_forward_alt < minimum_ratio:
                return False

            # reverse ratio test for the alternate allele
            if ratio_reverse_alt < minimum_ratio:
                return False

        return True
