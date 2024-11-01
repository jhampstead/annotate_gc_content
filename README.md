# annotate_gc_content
Annotates percentage GC content around VCF records.

## Usage
annotate_gc_content expects the following input arguments:
* A window size, describing the number of bases on either side of the variant to calculate percentage GC content
* A reference FASTA file
* An input VCF file
* An output VCF file

The program annotates an INFO/GC float tag in the output VCF file describing the proportion of GC bases within the specified window size (range 0 - 1).

## Compilation
HTSLib libraries are required to compile and execute. First load HTSLib, and then compile the executable using the following command:

```gcc -O3 -o annotate_mnv vcf_annotate_mnv.c -lhts```
