#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

// Function to calculate GC content around sequence variant
float calculate_gc_content(const char *sequence, int length) {
    int gc_count = 0;
    for (int i = 0; i < length; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C' || sequence[i] == 'g' || sequence[i] == 'c') {
            gc_count++;
        }
    }

    return (float)gc_count / length;
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <window_size> <reference.fasta> <input.vcf> <output.vcf>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int window_size = atoi(argv[1]);
    const char *fasta_file = argv[2];
    const char *vcf_file = argv[3];
    const char *output_file = argv[4];

    if (window_size <= 0) {
        fprintf(stderr, "Invalid window size: %d\n", window_size);
        return EXIT_FAILURE;
    }

    htsFile *in_vcf = bcf_open(vcf_file, "r");
    if (!in_vcf) {
        fprintf(stderr, "Could not open VCF file %s\n", vcf_file);
        return EXIT_FAILURE;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(in_vcf);
    if (!hdr) {
        fprintf(stderr, "Could not read VCF header\n");
        bcf_close(in_vcf);
        return EXIT_FAILURE;
    }

    // Add new INFO field for GC content
    int buffer_size = snprintf(NULL, 0, "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC content in %d base pair window around the variant\">\n", window_size) + 1; // +1 for the null terminator
    char *hdr_str = (char *)malloc(buffer_size);

    sprintf(hdr_str, "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC content in %d base pair window around the variant\">\n", window_size);
    bcf_hdr_append(hdr, hdr_str);
    htsFile *out_vcf = bcf_open(output_file, "w");
    bcf_hdr_write(out_vcf, hdr);

    free(hdr_str);

    // Open FASTA file
    faidx_t *fai = fai_load(fasta_file);
    if (!fai) {
        fprintf(stderr, "Could not open FASTA file %s\n", fasta_file);
        bcf_hdr_destroy(hdr);
        bcf_close(in_vcf);
        return EXIT_FAILURE;
    }

    // Read VCF records
    bcf1_t *rec = bcf_init();
    while (bcf_read(in_vcf, hdr, rec) == 0) {
        int pos = rec->pos; // 0-based position
        const char *chrom = bcf_seqname(hdr, rec);

        // Determine the window around the variant
        int start = pos - window_size / 2;
        int end = pos + window_size / 2;
        if (start < 0) start = 0;

        // Fetch the reference sequence
        int seq_len;
        char *seq = faidx_fetch_seq(fai, chrom, start, end, &seq_len);

        if (!seq) {
            fprintf(stderr, "Could not fetch sequence for %s:%d-%d\n", chrom, start, end);
            continue;
        }

        // Calculate the GC content
        float gc_content = calculate_gc_content(seq, seq_len);

        // Annotate VCF record
        bcf_update_info_float(hdr, rec, "GC", &gc_content, 1);
        bcf_write(out_vcf, hdr, rec);

        free(seq);
    }

    bcf_destroy(rec);
    fai_destroy(fai);
    bcf_hdr_destroy(hdr);
    bcf_close(in_vcf);
    bcf_close(out_vcf);

    return EXIT_SUCCESS;
}

