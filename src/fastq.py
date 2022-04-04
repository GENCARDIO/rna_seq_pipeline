#!/usr/bin/env python3

import sys
import re
import os.path
import logging
import pyfastx
from src.config import NGSbinaries

logger = logging.getLogger(__name__)

class FastqNotFound(Exception):
    pass

class MissingFastqPair(Exception):
    pass

class InvalidFastqNomenclature(Exception):
    '''
        Custom exception to track FASTQ invalid nomenclature
        :param str fq: fastq file
    '''

    def __init__(self, fq) -> None:
        self.message = "Invalid Fastq nomenclature: {}".format(fq)
        super().__init__(self.message)

class InvalidFastqFile(Exception):
    '''
        Custom exception to track FASTQ invalid files
        :param str fq: fastq file
    '''

    def __init__(self, fq) -> None:
        self.message = "Invalid Fastq file: {}".format(fq)
        super().__init__(self.message)

class Fastq:
    '''
        Fastq class. This class automatically validates fastq nomenclature and file consistency

        :param str fq: fastq file
        :py:meth:`check_consistency`
        :py:meth:`mean_read_length`
        :py:meth:`check_nomenclature`
    '''

    def __init__(self, fq, expect_paired=True, force_naming_convention=True):
        self._fq  = fq
        self._expect_paired = expect_paired
        self._paired = False

        if not os.path.isfile(self._fq):
            msg = ("Input Fastq {} not found").format(self._fq)
            logging.error(msg)
            raise FastqNotFound(msg)

        if self._expect_paired is True:
            if self.fq1 is not None and self.fq2 is not None:
                # print(self.fq1 + " " +self.fq2)
                self._paired = True
            else:
                msg = (" ERROR: Missing Fastq pair for {}").format(self._fq)
                raise MissingFastqPair(msg)

        if force_naming_convention is True:
            if self.check_nomenclature() is False:
                msg = (" ERROR: Invalid Fastq nomenclature for {}").format(self._fq)
                logging.error(msg)
                raise InvalidFastqNomenclature(self._fq)

    @property
    def fq1(self) -> str:
        '''
            :getter: Returns fastq1
        '''

        if '_R1_' in self._fq:
            return self._fq
        else:
            fq2 = self._fq
            fq1 = fq2.replace("R2", "R1")
            if self._expect_paired is True:
                if fq1 == fq2:
                    msg = (" ERROR: missing fastq pair for sample {}").format(self.sample_name)
                    logging.error(msg)
                    raise MissingFastqPair(msg)
        return(fq1)

    @property
    def fq2(self) -> str:
        '''
            :getter: Returns fastq2
        '''
        if '_R2_' in self._fq:
           return self._fq
        else :
           fq1 = self._fq
           fq2 = fq1.replace("R1", "R2")

           if self._expect_paired is True:
               if fq1 == fq2:
                 msg = (" ERROR: missing fastq pair for sample {}").format(self.sample_name)
                 logging.error(msg)
                 raise MissingFastqPair(msg)
        return(fq2)

    @property
    def fq1_basename(self) -> str:
        '''
            :getter: Returns fastq1 basename or prefix
        '''
        if self.fq1:
            fq1_basename = os.path.basename(self.fq1)
            return fq1_basename
        else:
            msg = (" ERROR: missing fastq1 {}").format(self.fq1)
            logging.error(msg)
            raise FileNotFoundError(msg)

    @property
    def fq2_basename(self) -> str:
        '''
            :getter: Returns fastq2 basename or prefix
        '''
        if self.fq2:
            fq2_basename = os.path.basename(self.fq2)
            return fq2_basename
        else:
            msg = (" ERROR: missing fastq2 {}").format(self.fq2)
            logging.error(msg)
            raise FileNotFoundError(msg)

    @property
    def sample_name(self) -> str:
        '''
            :getter: sample name from fastq prefix
            :var str fq: input fastq file
            :raises ValueError: if fastq name cannot be splitted by "_"
            :rtype: str
        '''

        fq_tmp = os.path.basename(self._fq).split("_")
        try:
            len(fq_tmp) > 1
        except:
            msg = (" ERROR: could not return sample name")
            logging.error(msg)
            raise ValueError()
        else:
            return(fq_tmp[0])

    @property
    def mean_read_length(self) -> int:
        '''
            Calculate the mean read length in base pairs

            :param fastq: fastq file to be checked
            :type: str
            :returns: returns the mean read length
            :rtype: int
        '''
        fq = pyfastx.Fastx(self._fq1)
        mean_length = fq.mean

        return mean_length

    def check_consistency(self) -> bool:
        '''
            Check that fastq1 and fastq2 are consistent:

            - Check that SEQ and QUAL have equal lengths
            - Check that fastq1 and fastq2 have equal read number

            :param str fastq1: fastq1 file to be checked
            :param str fastq2: fastq2 file to be checked

            :returns: returns True if no issues were detected, otherwise False
            :rtype: bool
        '''
        is_consistent = True
        fq1_reads = 0
        fq1 = pyfastx.Fastx(self.fq1)
        for name,seq,qual,comment in fq1:
            fq1_reads+=1
            if len(seq) != len(qual):
                msg = (" ERROR: Inconsistent length between SEQ and QUAL on line {} from file {}").\
                    format(str(fq1_reads), self.fq1)
                logging.error(msg)
                is_consistent = False

        if self._paired is True:
            fq2_reads = 0
            fq2 = pyfastx.Fastx(self.fq2)
            for name,seq,qual,comment in fq2:
                fq2_reads+=1
                if len(seq) != len(qual):
                    msg = (" ERROR: Inconsistent length between SEQ and QUAL on line {} from file {}").\
                        format(str(fq2_reads), self.fq2)
                    logging.error(msg)
                    is_consistent =  False

            if fq1_reads != fq2_reads:
                msg = (" ERROR: Unequal total reads between fq1 {}:{} and fq2 {}:{}").\
                    format( str(fq1_reads), self.fq1, str(fq2_reads), self.fq2)
                logging.error(msg)
                is_consistent = False
        return is_consistent

    def check_nomenclature(self) -> bool:
        '''
            Validate a fastq file:

            - Check FASTQ valid file extension (fq, fasta, fq.gz, fasta.gz)
            - Check that names are equal between fq1 and fq2,
            - Check illumina's nomenclature (_S[0-9]+_L[0-9]+_R[12]_[0-9]+)

            If paired:

            :param self.fq1: fastq1 to be checked
            :type: str
            :param self.fq2: fastq2 to be checked
            :type: str

            If single-ended:

            :param self.fq: fastq to be checked
            :type: str

            :returns: True if fastq(s) follow the Illumina's nomenclature
            :rtype: bool
        '''
        is_okay = True

        # For paired-end mode
        if self._paired is True:

            # First, Use regex to validate suffix
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self.fq1):
                pass
            else:
                is_okay = False
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self.fq2):
                pass
            else:
                is_okay = False

            # Second validate same names between fq1 and fq2
            fq1_tmp = os.path.basename(self.fq1).split("_")
            fq2_tmp = os.path.basename(self.fq2).split("_")

            fq1_name = fq1_tmp[0]
            fq2_name = fq2_tmp[0]

            if fq1_name == fq2_name:
                pass
            else:
                msg = (" ERROR: Inconsistent sample name between fastq1 {} and fastq2 {}").\
                    format(self.fq1, self.fq2)
                logging.error(msg)
                is_okay = False

            # Third, checking illumina's nomenclature
            if not re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self.fq1):
                msg = (" ERROR: Inconsistent Illumina's nomenclature on fastq1 {}").\
                    format(self.fq1)
                logging.error(msg)
                is_okay = False

            if not re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self.fq2):
                msg = (" ERROR: Inconsistent Illumina's nomenclature on fastq1 {}").\
                    format(self.fq2)
                logging.error(msg)
                is_okay = False
        # For single-ende mode
        elif self._paired is False:
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self._fq) is None:
                is_okay = False
            if re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self._fq) is None:
                is_okay = False
        return is_okay
