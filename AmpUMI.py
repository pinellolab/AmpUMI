import argparse
import sys
import os
import subprocess as sb
from sympy import *
from mpmath import *
import numpy as np
import re

class UMIException(Exception):
    pass

#def dedupUMIs(fastq,fastq_out,head_to_trim_seq,tail_to_trim_seq):
def dedupUMIs(args,parser):
    if not os.path.isfile(args.fastq):
        raise UMIException("Cannot find fastq file " + args.fastq)

    umi_regex_string = args.umi_regex

    #todo: collapse neighboring letters for better performance?
    umi_regex_string = umi_regex_string.replace('I','([ATCG])')
    umi_regex_string = umi_regex_string.replace('N','([ATCG])')
    umi_regex_string = umi_regex_string.replace('R','([AG])')
    umi_regex_string = umi_regex_string.replace('Y','([CT])')
    umi_regex_string = umi_regex_string.replace('S','([GC])')
    umi_regex_string = umi_regex_string.replace('W','([AT])')
    umi_regex_string = umi_regex_string.replace('K','([GT])')
    umi_regex_string = umi_regex_string.replace('M','([AC])')
    umi_regex_string = umi_regex_string.replace('B','([CGT])')
    umi_regex_string = umi_regex_string.replace('D','([AGT])')
    umi_regex_string = umi_regex_string.replace('H','([ACT])')
    umi_regex_string = umi_regex_string.replace('V','([ACG])')

    umi_regex_string = "(" + umi_regex_string + ")"

    umi_regex = re.compile(umi_regex_string)

    read_count = 0
    printed_count = 0
    count_with_regex = 0

    umi_keys_with_most_counts = {} #umi->key (key has most counts)
    umi_key_counts = {} #umi->count (count of fastqs key has seen) 

    umi_seq_counts = {} #umi_seq->count of that pair
    umi_seq_best_qual_fastqs = {} #umi_seq->fastq to be printed
    umi_seq_best_qual_sum = {} #best qual sum for the best fastq
    with open (args.fastq,"r") as f_in:
        while (1):
            id_line = f_in.readline()
            seq_line = f_in.readline()
            plus_line = f_in.readline()
            qual_line = f_in.readline()
            if not qual_line : break
            if not plus_line.startswith("+"):
                raise UMIException("Fastq %s cannot be parsed (%s%s%s%s) "%(args.fastq,id_line,seq_line,plus_line,qual_line))
            read_count += 1
            #print('seq:'+seq_line)
            #print('umi: ' + umi_regex_string)

            match_obj = umi_regex.search(seq_line.upper())
            if match_obj:
                count_with_regex += 1

                #group 1 is the whole string
                this_UMI = ""
                for i in range(2,len(match_obj.groups())+1):
                    #print("got group " + str(i) + ": " + match_obj.group(i))
                    #print("start " + str(match_obj.start(i)))
                    #print("end " + str(match_obj.end(i)))
                    this_UMI = this_UMI + match_obj.group(i)

                trimmed_seq = seq_line[0:match_obj.start(1)] + seq_line[match_obj.end(1):]
                trimmed_qual = qual_line[0:match_obj.start(1)] + qual_line[match_obj.end(1):]
                trimmed_qual_sum = np.sum(np.fromstring(trimmed_qual,dtype=np.uint8))

                this_key = this_UMI + " # " + trimmed_seq


                if this_key not in umi_seq_counts:
                    umi_seq_counts[this_key] = 1
                    umi_seq_best_qual_sum[this_key] = trimmed_qual_sum
                    umi_seq_best_qual_fastqs[this_key] = id_line + trimmed_seq + plus_line + trimmed_qual
                else:
                    umi_seq_counts[this_key] += 1
                    if umi_seq_best_qual_sum[this_key] < trimmed_qual_sum:
                        umi_seq_best_qual_sum[this_key] = trimmed_qual_sum
                        umi_seq_best_qual_fastqs[this_key] = id_line + trimmed_seq + plus_line + trimmed_qual

                if this_UMI not in umi_key_counts:
                    umi_key_counts[this_UMI] = 1
                    umi_keys_with_most_counts[this_UMI] = this_key
                else:
                    if umi_seq_counts[this_key] > umi_key_counts[this_UMI]:
                        umi_key_counts[this_UMI] =  umi_key_counts[this_UMI]
                        umi_keys_with_most_counts[this_UMI] = this_key

    if read_count == 0:
        raise UMIException("UMI command failed. Got no reads from " + args.fastq )

    print("Read " + str(read_count) + " reads from " + args.fastq)
    print("Of those, " + str(count_with_regex) + " had the regex")

    umi_list = sorted(umi_key_counts, key=lambda k: umi_key_counts[k])
    umi_count = len(umi_list)
#    print("umi_list: " + str(umi_list))

    print("Processed " + str(umi_count) + " barcodes")
    if umi_count == 1 and umi_list[0] == '':
        print("Warning, only the empty barcode '' was found.")

    f_out = open(args.fastq_out,"w")


    collision_count = 0
    collision_count_reads = 0
    too_few_reads_count = 0
    too_few_reads_count_reads = 0
    printed_count = 0
    printed_count_reads = 0
    for umi_seq in umi_seq_counts:
        umi, seq = umi_seq.split(" # ")
        #if another umi had more reads than ths one...
        if umi_keys_with_most_counts[umi] != umi_seq:
            collision_count += 1
            collision_count_reads += umi_seq_counts[umi_seq]
        #if this umi had too few reads
        elif umi_seq_counts[umi_seq] < args.min_UMI_to_print:
            too_few_reads_count += 1
            too_few_reads_count_reads += umi_seq_counts[umi_seq]
        else:
            f_out.write(umi_seq_best_qual_fastqs[umi_seq])
            printed_count += 1
            printed_count_reads += umi_seq_counts[umi_seq]

    f_out.close()

    print("Observed %d UMI collisions (%d reads). Not printing these."%(collision_count,collision_count_reads))
    print("Observed %d UMIs (%d reads) with too few reads. Not printing these."%(too_few_reads_count,too_few_reads_count_reads))
    print("Printed %d deduplicated sequences (%d reads) to %s"%(printed_count,printed_count_reads,args.fastq_out))


def calculateUMIsmath(numUMIs,numMolecules):
        try:
                xrange
        except NameError:
                xrange = range

        p=1.0
        for i in xrange(1,numMolecules):
                p*=(1.0-i/numUMIs)
        return p

def calculateUMIs(args,parser):
            umiCount = 0
            umiLength = 0
            if args.nu:
                umiCount = args.nu
                umiLength = mpmath.log(umiCount,4)
                Q = mpf(umiCount)
            elif args.ul:
                umiCount = 4**args.ul
                umiLength = args.ul
                Q = mpf(umiCount)
            else:
                parser.print_help()
                exit("Either -ul or -nu is required for umi length calculation")
            p = calculateUMIsmath(Q,args.nm)
            print("With %d umis (length %d) and %d unique molecules, the probability of no collisions is %f"%(umiCount,umiLength,args.nm,p))

def main():
        parser = argparse.ArgumentParser(prog='AmpUMI - A toolkit for designing and analyzing amplicon sequencing experiments using unique molecular identifiers')

        subparsers = parser.add_subparsers(help='Enter a specific AmpUMI function',dest='subparser_name')

        parser_run = subparsers.add_parser('Process',help='Process a fastq with UMIs for downstream processing')
        parser_run.add_argument('--fastq',required=True,help="Path to the fastq to be processed")
        parser_run.add_argument('--fastq_out',required=True,help="Path to the trimmed fastq to be written")
        parser_run.add_argument('--umi_regex',help='Regular expression specifying the umi (I) as well as any primer sequences to be trimmed (A,C,T,G)\nFor example, if the UMI is the first 5 basepairs, this should be "^IIIII".',required=True)
        parser_run.add_argument('--min_UMI_to_print',help='The minimum times a UMI must be seen to be printed',type=int,default=0)
        parser_run.set_defaults(func=dedupUMIs)

        parser_calculate = subparsers.add_parser('Calculate',help='Calculate UMI collision probability')

        parser_calculate.add_argument("-ul",type=int,help="UMI length",required=false)
        parser_calculate.add_argument("-nu",type=int,help="Number of unique umis",required=false)
        parser_calculate.add_argument("-nm",type=int,help="Number of unique molecules")
        parser_calculate.set_defaults(func=calculateUMIs)

        if len(sys.argv)==1:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args()
        args.func(args,parser)

#        dedupUMIs(args.fastq,args.fastq_out,args.head_to_trim_seq,args.tail_to_trim_seq)


if __name__ == "__main__":
    main()
