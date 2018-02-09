# AmpUMI
Toolkit for the design and analysis of amplicon sequencing experiments utilizing unique molecular identifiers (UMIs).

## Usage:
After downloading AmpUMI.py, the program can be run to and command line parameters can be shown using the command:
```
python AmpUMI.py -h
```
AmpUMI can be run in two modes: *Calculate* and *Process*.
```
python AmpUMI.py Calculate -h
python AmpUMI.py Process -h
```

### Calculate mode
The calculate mode is used to compute the probability of a UMI/molecule collision given a number of unique UMIs or the UMI length and given a number of molecules that will be paired with the UMIs. We calculate this probability for the worst case in which all molecules have the same sequence.

To calculate this probability, AmpUMI Calculate should be run with the following parameters:
*  ```-ul``` UMI length
*  ```-nu``` Number of unique UMIs
*  ```-nm``` Number of molecules

Either ```-nu``` or ```-ul``` should be provided. If ```-ul``` is provided, ```-nu``` will be calculated as ```4^ul```.

For example, if we are using an 15bp UMI, and plan to add 10,000 molecules to pair with the UMIs, the probability of observing a collision can be calculated with the following command:
```
python AmpUMI.py Calculate -ul 15 -nm 10000
```
This produces the following output:
```
With 1073741824 umis (length 15) and 10000 unique molecules, the probability of no collisions is 0.954506
```

### Process mode
The process mode is used to process amplicon sequencing libraries containing UMIs. 

AmpUMI Process is run with the following parameters:
*  ```--fastq``` Path to the fastq to be processed
*  ```--fastq_out``` Path to the trimmed fastq to be written
*  ```--umi_regex``` Regular expression specifying the umi (I) as well as any primer sequences to be trimmed (A,C,T,G).
*  ```--min_UMI_to_print``` The minimum times a UMI must be seen to be printed (default=0)

AmpUMI Process mode will parse and trim the UMi and any specified adapter from each read. Next, for each UMI-sequence pair, the highest-quality sequence will be kept. Next, error correction will be performed in two steps: 1) for each UMI, the most prevalent sequence will be kept, and other less-frequent sequences will be discarded. Thus, if there are sequencing errors in the read, these will be filtered out. 2) each UMI is printed only if it was seen at least ```--min_UMI_to_print``` times. This is good for filtering sequencing errors that may affect the UMI sequence. In this case, UMIs with few reads may be the result of sequencing errors in those UMIs. Finally, duplicate reads are removed. For each UMI, only one read is printed to the final output, removing PCR duplicates. 

For example, if the UMI is the first 5 basepairs of a read, and the UMI consists of any possible base combination, the following command should be used:
```
python AmpUMI.py --fastq input.fastq --fastq_out input.fastq.dedup.fastq --umi_regex "^IIIII"
```
In this example, the *^* symbol anchors the UMI match to the beginning of the read, and the next 5bp are the UMI. 

If the last 5bp of each read were the UMI, this could be specified using the flag
```--umi_regex "IIIII$```
where the *$* symbol anchors the UMI match to the beginning of the read. 

[Ambiguous characters](https://www.bioinformatics.org/sms/iupac.html) can be used to specify UMI locations as well. 

If adapter sequences are present in reads, or surround UMIs, AmpUMI will trim these adapter sequences as well. For example, if your UMI design is a 5bp UMI flanked by the adapter sequnce ACCTG, this could be specified using the flag
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
The ```--min_UMI_to_print``` parameter sets the requirement for the minimum number of reads that must be seen for a particular UMI for it to be printed.




