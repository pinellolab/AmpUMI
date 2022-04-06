# AmpUMI
Toolkit for the design and analysis of amplicon sequencing experiments utilizing unique molecular identifiers (UMIs). 

For details, please see our publication at https://doi.org/10.1093/bioinformatics/bty264.

## Installation:
AmpUMI requires the sympy, mpmath, and numpy packages, and can be installed using the command:
```
pip install git+https://github.com/pinellolab/AmpUMI.git
```

## Usage:
After downloading AmpUMI.py, the program can be run and command line parameters can be shown using the command:
```
python AmpUMI.py -h
```
AmpUMI can be run in four modes: *Collision*, *CollisionNumber*, *Distortion* and *Process*. The *Collision*, *CollisionNumber* and *Distortion* commands are for use in selecting an appropriate UMI length during the experiment design process, and *Process* is used for processing FASTQ reads after an amplicon sequencing experiment.
```
python AmpUMI.py Collision -h
python AmpUMI.py CollisionNumber -h
python AmpUMI.py Distortion -h
python AmpUMI.py Process -h
```

### Collision mode
The *Collision* mode is used to compute the probability of having no UMI/molecule collisions given a number of unique UMIs or the UMI length and given a number of molecules that will be paired with the UMIs. We calculate this probability for the worst case in which all molecules have the same sequence.

To calculate this probability, AmpUMI Collision should be run with the following parameters:
*  ```-ul``` UMI length or ```-nu``` Number of unique UMIs
*  ```-nm``` Number of molecules

Either ```-nu``` or ```-ul``` should be provided. If ```-ul``` is provided, ```-nu``` will be calculated as ```4^ul```.

For example, if we are using a 15bp UMI, and plan to add 10,000 molecules to pair with the UMIs, the probability of observing no collisions can be calculated with the following command:
```
python AmpUMI.py Collision -ul 15 -nm 10000
```
This produces the following output:
```
With 1073741824 UMIs (length 15) and 10000 unique molecules, the probability of no collisions is 0.954506
```

The minimum barcode length required to have a probability of collisions can also be calculated using this command by using the ```--min_collision_p``` parameter. For example, to compute the minimum UMI length in order to observe no collisions in an experiment with 10,000 molecules, AmpUMI should be run as follows:
```
python AmpUMI.py Collision -nm 10000 -mp 0.95
```
This produces the following output:
```
With 1073741824 UMIs (length 15) and 10000 unique molecules, the probability of no collisions is 0.954506
```

The minimum barcode length required to have a probability of collisions can also be calculated using this command by using the ```--min_collision_p``` parameter. For example, to compute the minimum UMI length in order to observe no collisions in an experiment with 10,000 molecules, AmpUMI should be run as follows:
```
python AmpUMI.py Collision -nm 10000 -mp 0.95
```
This produces the following output:
```
With 1073741824 UMIs (length 15) and 10000 unique molecules, the probability of no collisions is 0.954506
```



### CollisionNumber mode
The *CollisionNumber* mode is used to calculate the expected number of collisions given a UMI length and the number of molecules in the experiment, and the distribution of alleles in the experiment.

To calculate the number of collisions, AmpUMI CollisionNumber should be run with the following parameters:
*  ```-ul``` UMI length or ```-nu``` Number of unique UMIs
*  ```-nm``` Number of molecules
*  ```-af``` Comma-separated list of allele frequencies or allele fractions

Either ```-nu``` or ```-ul``` should be provided. If ```-ul``` is provided, ```-nu``` will be calculated as ```4^ul```.

The case with the highest expected number of collisions is the case where there is only one allele. If the allele distribution is not known, this worst case should be calculated by running with the argument ```-af 1```.  

For example, if we are using a 15bp UMI, and plan to add 10,000 molecules to pair with the UMIs, the expected number of collisions can be calculated with the following command:

```
python AmpUMI.py CollisionNumber -ul 15 -nm 10000 -af 1
```

This produces the following output:
```
With 1073741824 UMIs (length 15) and 10000 molecules, the expected number of collisions is 0.046561

Allelic fraction  Number of collisions
1                 0.0465613603591919

```

Thus, there is fewer than one (0.047) collision expected for this experimental setup. For each allele fraction specified with the ```-af``` flag, the output describes how many collisions are expected for that allele. In this case, there is only one allele (100% of alleles) and there are 0.047 collisions expected for that allele. 

With more than one allele, if we were using a 8bp UMI with 10,000 molecules, with a 40/60 mixture of alleles, we could calculate the expected number of collisions using the command: 

```
python AmpUMI.py CollisionNumber -ul 8 -nm 10000 -af 0.4,0.6
```

This produces the following output:
```
With 65536 UMIs (length 8) and 10000 molecules, the expected number of collisions is 386.052363

Allelic fraction      Number of collisions
0.4                   119.612739004733
0.6                   266.439623527411

```

The minimum barcode length required to have fewer than a specified number of expected collisions can also be calculated using this command by using the ```--max_collision_number``` (```-mn```)  parameter. For example, to compute the minimum UMI length in order to observe no greater than 5 expected collisions in an experiment with 10,000 molecules, AmpUMI should be run as follows:
```
python AmpUMI.py CollisionNumber -nm 10000 -mn 5 -af 1
```
This produces the following output:
```
With 16777216 UMIs (length 12) and 10000 molecules, the expected number of collisions is 2.979342

Allelic fraction       Number of collisions
1                      2.97934236191213
```


### Distortion mode
In amplicon sequencing experiments using UMIs, the observed allelic fraction of each variant can be distorted if the UMI length is too short and the UMI complexity is too low. In the case in which allele variants are present at roughly the same proportion in the population this distortion will be small, as UMI-molecule collisions will affect each molecule type equally. However, as the proportion of allele variants becomes more imbalanced, after deduplication using UMIs the observed frequency of rare alleles will increase, and the observed frequency of frequent alleles will be lower than the actual frequency in the population. Intuitively, this is because collisions are more frequent in the more frequent alleles, where a larger UMI diversity is required to avoid UMI-molecule collisions. Rare allele with few reads are less likely to exhaust the available unique UMIs. The *Distortion* mode can be used to compute the expected allelic distortion that would arise from having a specific UMI length.


To calculate this expected distortion, AmpUMI Distortion should be run with the following parameters:
*  ```-af``` Comma-separated list of allele frequencies or allele fractions
*  ```-ul``` UMI length or ```-nu``` Number of unique UMIs
*  ```-nm``` Number of molecules

Either ```-nu``` or ```-ul``` should be provided. If ```-ul``` is provided, ```-nu``` will be calculated as ```4^ul```.

For example, in an experiment where the allele frequencies are 50%, 30%, 10%, and 10%, with a UMI length of 5 and 1000 molecules, the expected allelic distortion can be computed using the command:
```
python AmpUMI.py Distortion -af 0.5,0.3,0.1,0.1 -nm 1000 -ul 5
```
This produces the following output:
```
With 1024 UMIs (length 5) and 1000 molecules, the expected total allelic fraction distortion is 0.064952
Actual  Expected after deduplication
0.5     0.467523770626148
0.3     0.307316077591381
0.1     0.112580075891235
0.1     0.112580075891235
```
This demonstrates that the presence of minor alleles (i.e. the alleles with 10% presence in the real population) are overestimated with barcode complexity that is too low, and that major alleles (i.e. the allele with 50% presence in the real population) are underestimated. 

The minimum barcode length to have an allelic fraction distortion below a certain threshold can be computed by using the ```--max_distortion``` parameter. For example, in the experiment above, with allele frequencies of 50%, 30%, 10%, and 10%, the minimum barcode length to have a total allelic distortion below 1% would be computed using the following command:
```
python AmpUMI.py Distortion -af 0.5,0.3,0.1,0.1 -nm 1000 --max_distortion 0.01
```
This produces the following output:
```
With 16384 UMIs (length 7) and 1000 molecules, the expected total allelic fraction distortion is 0.004256
Actual   Expected after deduplication
0.5      0.497872213416257
0.3      0.300542912015644
0.1      0.10079243728405
0.1      0.10079243728405
```

### Process mode
The process mode is used to process amplicon sequencing libraries containing UMIs. 

AmpUMI Process is run with the following parameters:
*  ```--fastq``` Path to the fastq to be processed
*  ```--fastq_out``` Path to the trimmed fastq to be written
*  ```--umi_regex``` Regular expression specifying the umi (I) as well as any primer sequences to be trimmed (A,C,T,G).
*  ```--min_umi_to_keep``` The minimum times a UMI must be seen to be keep (default=0)
*  ```--write_UMI_counts``` Flag to write counts of each UMI to a file
*  ```--write_alleles_with_multiple_UMIs``` Flag to write alleles with multiple UMIs to a file


AmpUMI Process mode will parse and trim the UMI and any specified adapter from each read. Next, for each UMI-sequence pair, the highest-quality sequence will be kept. Next, error correction will be performed in two steps: 1) for each UMI, the most prevalent sequence will be kept, and other less-frequent sequences will be discarded. Thus, if there are sequencing errors in the read, these will be filtered out. 2) (optional) each UMI is printed only if it was seen at least ```--min_umi_to_print``` times. This is good for filtering sequencing errors that may affect the UMI sequence. In this case, UMIs with few reads may be the result of sequencing errors in those UMIs. Finally, duplicate reads are removed. For each UMI, only one read is printed to the final output, removing PCR duplicates. 

For example, if the UMI is the first 5 basepairs of a read, and the UMI consists of any possible base combination, the following command should be used:
```
python AmpUMI.py --fastq input.fastq --fastq_out input.fastq.dedup.fastq --umi_regex "^IIIII"
```
In this example, the *^* symbol anchors the UMI match to the beginning of the read, and the next 5bp are the UMI. 

If the last 5bp of each read were the UMI, this could be specified using the flag
```--umi_regex "IIIII$```
where the *$* symbol anchors the UMI match to the end of the read. 

[Ambiguous characters](https://www.bioinformatics.org/sms/iupac.html) can be used to specify UMI locations as well. 

If adapter sequences are present in reads, or surround UMIs, AmpUMI will trim these adapter sequences as well. For example, if your UMI design is a 5bp UMI flanked by the adapter sequence ACCTG, this could be specified using the flag
```--umi_regex "ACCTGIIIIIACCTG"```

In addition, if the adapter or UMI are not always exactly at the beginning of the read, bases before the adapter and UMI (as well as the adapter and UMI) will be trimmed using by using the .* operator which represents any character. 

For example:
```
--umi_regex ".*ACCTGIIIIIACCTG"
pre-AmpUMI:  TTTACCTGAAAAAACCTGTATATATAT
post-AmpUMI:                   TATATATAT (UMI: AAAAA)
pre-AmpUMI:  TACCTGATATAACCTGGAGAGAG
post-AmpUMI:                 GAGAGAG (UMI: ATATA)
```
The ```--min_umi_to_print``` parameter sets the requirement for the minimum number of reads that must be seen for a particular UMI for it to be printed. Setting this parameter will discard UMI-read pairs that occur less than a specified number of times. This can remove spurious UMIs or reads that arise from sequencing error and only appear infrequently as opposed to the more common variants that arise from source biological material.

Setting the ```--write_UMI_counts``` flag writes the number of times each UMI appeared in the sample. If a few barcodes dominate the pool, there may have been a problem with the original UMI pool or with the library construction or amplification process. 

Setting the ```--write_alleles_with_multiple_UMIs``` flag writes UMI/allele pairs for which the UMI was paired to alleles with different sequences. These are UMI collisions and these reads are discarded from the main output. 

## Citation
If you use AmpUMI, please cite our publication: 

Kendell Clement, Rick Farouni, Daniel E Bauer, Luca Pinello; AmpUMI: design and analysis of unique molecular identifiers for deep amplicon sequencing, Bioinformatics, Volume 34, Issue 13, 1 July 2018, Pages i202–i210, https://doi.org/10.1093/bioinformatics/bty264


