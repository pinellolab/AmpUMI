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
                    #if this sequence has the highest quality, store it
                    if umi_seq_best_qual_sum[this_key] < trimmed_qual_sum:
                        umi_seq_best_qual_sum[this_key] = trimmed_qual_sum
                        umi_seq_best_qual_fastqs[this_key] = id_line + trimmed_seq + plus_line + trimmed_qual

                if this_UMI not in umi_key_counts:
                    umi_key_counts[this_UMI] = 1
                    umi_keys_with_most_counts[this_UMI] = this_key
                else:
                    umi_key_counts[this_UMI] += 1
                    #if this sequence is the most seen for this UMI, store it
                    if umi_seq_counts[this_key] > umi_key_counts[this_UMI]:
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
        raise Exception("Error: only the empty barcode '' was found.")

    f_out = open(args.fastq_out,"w")
    if (args.write_alleles_with_multiple_UMIs):
        f_out_multiUMI = open(args.fastq_out+".AmpUMI.multi.out","w")
        f_out_multiUMI.write('discarded/kept\tCount\tUMI\tSequence\n')

    collision_count = 0
    collision_count_reads = 0
    too_few_reads_count = 0
    too_few_reads_count_reads = 0
    printed_count = 0
    printed_count_reads = 0
    collided_umi_reads = []
    for umi_seq in umi_seq_counts:
        umi, seq = umi_seq.split(" # ")
        #if another umi had more reads than ths one...
        if umi_keys_with_most_counts[umi] != umi_seq:
            collision_count += 1
            collision_count_reads += umi_seq_counts[umi_seq]
            if (args.write_alleles_with_multiple_UMIs):
                other_key = umi_keys_with_most_counts[umi]
                other_umi, other_seq = other_key.split(" # ")
                f_out_multiUMI.write('k\t' + str(umi_seq_counts[other_key]) + "\t" + other_umi + "\t" + other_seq)
                f_out_multiUMI.write("d\t" + str(umi_seq_counts[umi_seq])+"\t" + umi + "\t" + seq)
        #if this umi had too few reads
        elif umi_seq_counts[umi_seq] < args.min_umi_to_keep:
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

    if (args.write_UMI_counts):
        f_out_ampUMI = open(args.fastq_out+".AmpUMI.out","w")
        f_out_ampUMI.write('UMI\tCount\n')
        for umi in sorted(umi_key_counts):
            f_out_ampUMI.write(umi + "\t" + str(umi_key_counts[umi]) + "\n")
        print('Wrote UMI counts to ' + args.fastq_out+'.AmpUMI.out')

    if (args.write_alleles_with_multiple_UMIs):
        print('Wrote UMI collisions to ' + args.fastq_out+'.AmpUMI.multi.out')


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
            if (args.nm is None):
                parser.print_help()
                exit("Molecule count -nm is required")
            # determine min umi len to have this max distortion
            if args.mp is not None:
                if args.mp == 1:
                    exit("Non-collision p-value of 1 cannot be calculated")
                for i in range(50):
                    umi_count = 4**i
                    p = calculateUMIsmath(mpf(umi_count),args.nm)
                    if p >= args.mp:
                        print("With %d UMIs (length %d) and %d unique molecules, the probability of no collisions is %f"%(umi_count,i,args.nm,p))
                        exit()
                raise Exception('Cannot find barcode length producing a p-value greater than %f', args.mp)

            umiCount = 0
            umiLength = 0
            if args.nu:
                umiCount = args.nu
                umiLength = log(umiCount,4)
                Q = mpf(umiCount)
            elif args.ul:
                umiCount = 4**args.ul
                umiLength = args.ul
                Q = mpf(umiCount)
            else:
                parser.print_help()
                exit("Either -ul or -nu is required for umi length calculation")

            p = calculateUMIsmath(Q,args.nm)
            print("With %d UMIs (length %d) and %d unique molecules, the probability of no collisions is %f"%(umiCount,umiLength,args.nm,p))

def calculateCollisionNumberMath(allele_fracs,umi_count,molecule_count):
        var_I = len(allele_fracs)
        var_J = mpf(umi_count)
        var_n = mpf(molecule_count)
        collision_counts = []
        sum_collisions = 0
        for i in range(len(allele_fracs)):
            observed_m = mpf(allele_fracs[i])
            this_c = var_J * (1-exp(var_n *log(1-(observed_m/var_J))))
            this_collisions = observed_m * var_n - this_c

            sum_collisions += this_collisions
            collision_counts.append(this_collisions)

        return sum_collisions,collision_counts

def calculateCollisionNumber(args,parser):
            if args.af is None:
                parser.print_help()
                exit("Allele frequencies or allele fractions -af are required")
            if args.nm is None:
                parser.print_help()
                exit("Molecule count -nm is required")

            parsed_allele_fracs = args.af.split(',')
            allele_fracs = parsed_allele_fracs
            is_allele_fractions = all(mpf(x) <=1 for x in parsed_allele_fracs)
            if not is_allele_fractions:
                allele_fracs = [mpf(x)/args.ns for x in parsed_allele_fracs]

            # determine min umi len to have this max number of collisions
            if args.mn is not None:
                for i in range(50):
                    umi_count = 4**i
                    sum_collisions,collision_counts = calculateCollisionNumberMath(allele_fracs = allele_fracs,umi_count =
                            umi_count,molecule_count = args.nm)
                    if sum_collisions <= args.mn:
                        print("With %d UMIs (length %d) and %d molecules, the expected number of collisions is %f"%(umi_count,i,args.nm,sum_collisions))
                        print("Allelic fraction\tNumber of collisions")
                        for i in range(len(allele_fracs)):
                            print(str(allele_fracs[i])+ '\t' + str(collision_counts[i]))
                        exit()
                raise Exception('Cannot find barcode length producing fewer than %f collisions', args.mp)

            #otherwise, calculate distortion
            umiCount = 0
            umiLength = 0
            if args.nu:
                umiCount = args.nu
                umiLength = log(umiCount,4)
                Q = mpf(umiCount)
            elif args.ul:
                umiCount = 4**args.ul
                umiLength = args.ul
                Q = mpf(umiCount)
            else:
                parser.print_help()
                exit("Either -ul or -nu is required")


            sum_collisions,collision_counts = calculateCollisionNumberMath(allele_fracs = allele_fracs,umi_count = Q,molecule_count = args.nm)
            print("With %d UMIs (length %d) and %d molecules, the expected number of collisions is %f"%(umiCount,umiLength,args.nm,sum_collisions))
            print("Allelic fraction\tNumber of collisions")
            for i in range(len(allele_fracs)):
                print(str(allele_fracs[i])+ '\t' + str(collision_counts[i]))

def calculateDistortionMath(allele_fracs,umi_count,molecule_count):
        var_I = len(allele_fracs)
        var_J = mpf(umi_count)
        var_n = mpf(molecule_count)
        predicted_m = []
        sum_predicted = 0
        for observed_m in allele_fracs:
            this_m = 1-exp(var_n *log(1-(mpf(observed_m)/var_J)))
            sum_predicted += this_m
            predicted_m.append(this_m)

        outer = 1/sum_predicted

        sum_allelic_distortion = 0

        final_allele_fracs = []
        for i in range(len(allele_fracs)):
            real_m = allele_fracs[i]
            pred_m = outer * predicted_m[i]
            final_allele_fracs.append(pred_m)
            sum_allelic_distortion += abs(mpf(real_m) - pred_m)
        return sum_allelic_distortion,final_allele_fracs

def calculateDistortion(args,parser):
            if args.af is None:
                parser.print_help()
                exit("Allele frequencies or allele fractions -af are required")
            if args.nm is None:
                parser.print_help()
                exit("Molecule count -nm is required")

            parsed_allele_fracs = args.af.split(',')
            allele_fracs = parsed_allele_fracs
            is_allele_fractions = all(mpf(x) <=1 for x in parsed_allele_fracs)
            if not is_allele_fractions:
                allele_fracs = [mpf(x)/args.ns for x in parsed_allele_fracs]

            # determine min umi len to have this max distortion
            if args.md is not None:
                if args.md == 0:
                    exit("Max distortion of 0 cannot be calculated")
                for i in range(50):
                    umi_count = 4**i
                    distortion,new_allele_fracs = calculateDistortionMath(allele_fracs = allele_fracs,umi_count =
                            umi_count,molecule_count = args.nm)
                    if distortion <= args.md:
                        print("With %d UMIs (length %d) and %d molecules, the expected total allelic fraction distortion is %f"%(umi_count,i,args.nm,distortion))
                        print("Actual\tExpected after deduplication")
                        for i in range(len(allele_fracs)):
                            print(str(allele_fracs[i])+ '\t' + str(new_allele_fracs[i]))
                        exit()
                raise Exception('Cannot find barcode length producing a distortion less than %f', args.mp)
                
            #otherwise, calculate distortion
            umiCount = 0
            umiLength = 0
            if args.nu:
                umiCount = args.nu
                umiLength = log(umiCount,4)
                Q = mpf(umiCount)
            elif args.ul:
                umiCount = 4**args.ul
                umiLength = args.ul
                Q = mpf(umiCount)
            else:
                parser.print_help()
                exit("Either -ul or -nu is required")


            distortion,new_allele_fracs = calculateDistortionMath(allele_fracs = allele_fracs,umi_count = Q,molecule_count = args.nm)
            print("With %d UMIs (length %d) and %d molecules, the expected total allelic fraction distortion is %f"%(umiCount,umiLength,args.nm,distortion))
            print("Actual\tExpected after deduplication")
            for i in range(len(allele_fracs)):
                print(str(allele_fracs[i])+ '\t' + str(new_allele_fracs[i]))

def main():
        #parser = argparse.ArgumentParser(prog='AmpUMI - A toolkit for designing and analyzing amplicon sequencing experiments using unique molecular identifiers\n')
        parser = argparse.ArgumentParser(description='AmpUMI - A toolkit for designing and analyzing amplicon sequencing experiments using unique molecular identifiers\n')
        parser.add_argument('--version', action='version', version='%(prog)s 1.2')

        subparsers = parser.add_subparsers(help='Enter a specific AmpUMI function',dest='subparser_name')

        #PROCESS
        parser_run = subparsers.add_parser('Process',help='Process a fastq with UMIs for downstream processing')
        parser_run.add_argument('--fastq',required=True,help="Path to the fastq to be processed")
        parser_run.add_argument('--fastq_out',required=True,help="Path to the trimmed fastq to be written")
        parser_run.add_argument('--umi_regex',help='Regular expression specifying the umi (I) as well as any primer sequences to be trimmed (A,C,T,G)\nFor example, if the UMI is the first 5 basepairs, this should be "^IIIII".',required=True)
        parser_run.add_argument('--min_umi_to_keep',help='The minimum times a UMI must be seen to be kept',type=int,default=0)
        parser_run.add_argument('--write_UMI_counts',help='Flag to write counts of each UMI to a file',action='store_true')
        parser_run.add_argument('--write_alleles_with_multiple_UMIs',help='Flag to write alleles with multiple UMIs to a file',action='store_true')
        parser_run.set_defaults(func=dedupUMIs)

        #CALCULATE COLLISION
        parser_collision = subparsers.add_parser('Collision',help='Calculate UMI collision probability')
        parser_collision_group = parser_collision.add_mutually_exclusive_group(required=True)
        parser_collision.add_argument("-nm","--number_molecules",dest="nm",type=int,help="Number of unique molecules",required=True)
        parser_collision_group.add_argument("-ul","--UMI_Length",dest="ul",type=int,help="UMI length",required=false)
        parser_collision_group.add_argument("-nu","--number_UMIs",dest="nu",type=int,help="Number of unique UMIs",required=false)
        parser_collision_group.add_argument("-mp","-min_collision_p",dest="mp",type=float,help="Minimum collision probability (If this argument is given, the program returns the minimum barcode length for which a probability of observing no collisions is at least this value.")
        parser_collision.set_defaults(func=calculateUMIs)

        #DISTORTION
        parser_distortion = subparsers.add_parser('Distortion',help='Calculate distortion of real vs observed allele frequencies')
        parser_distortion_group = parser_distortion.add_mutually_exclusive_group(required=True)
        parser_distortion.add_argument("-af",help='Comma-separated list of real allele frequencies or allele fractions',required=True)
        parser_distortion.add_argument("-nm","--number_molecules",dest="nm",type=int,help="Number of molecules",required=True)
        parser_distortion_group.add_argument("-ul","--UMI_Length",dest="ul",type=int,help='UMI length')
        parser_distortion_group.add_argument("-nu","--number_UMIs",dest="nu",type=int,help='Number of unique UMIs')
        parser_distortion_group.add_argument("-md","--max_distortion",dest="md",type=float,help="Maximum distortion (If this argument is given, the program returns the minimum barcode length for which the distortion is smaller than this value.")
        parser_distortion.set_defaults(func=calculateDistortion)

        #COLLISION NUMBER
        parser_collision_number = subparsers.add_parser('CollisionNumber',help='Calculate expected number of collisions')
        parser_collision_number_group = parser_collision_number.add_mutually_exclusive_group(required=True)
        parser_collision_number.add_argument("-af",help='Comma-separated list of real allele frequencies or allele fractions',required=True)
        parser_collision_number.add_argument("-nm","--number_molecules",dest="nm",type=int,help="Number of molecules",required=True)
        parser_collision_number_group.add_argument("-ul","--UMI_Length",dest="ul",type=int,help='UMI length')
        parser_collision_number_group.add_argument("-nu","--number_UMIs",dest="nu",type=int,help='Number of unique UMIs')
        parser_collision_number_group.add_argument("-mn","--max_collision_number",dest="mn",type=float,help="Maximum collision number (If this argument is given, the program returns the minimum barcode length for which the number of collisions is smaller than this value.")
        parser_collision_number.set_defaults(func=calculateCollisionNumber)

        if len(sys.argv)==1:
            parser.print_help()
            sys.exit(1)
        args = parser.parse_args()
        args.func(args,parser)

#        dedupUMIs(args.fastq,args.fastq_out,args.head_to_trim_seq,args.tail_to_trim_seq)


if __name__ == "__main__":
    main()
