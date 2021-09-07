# aliquot-maf-tools

Tools for creating, merging, and filtering aliquot-level GDC MAFs. These tools depend on the 
[maf-lib](https://github.com/NCI-GDC/maf-lib) library and requires strict schema definitions.
There are three main utilities in the `aliquot-maf-tools` software:

* `VcfToAliquotMaf` - Converts a single GDC annotated VCF to a caller-specific aliquot MAF
* `MergeAliquotMafs` - Takes 2 or more caller-specific aliquot MAFs for a single tumor/normal pair and merges into a merged multi-caller MAF
* `MaskMergedAliquotMaf` - Filters the merged multi-caller MAF

## Convert VCF to Aliquot MAF

All of the implementations use the format:

```
aliquot-maf-tools VcfToAliquotMaf \
    --input_vcf <path to input annotated VCF> \
    --output_maf <path to output MAF> \
    <schema version> <OPTIONS>
```

We expect the VCF file to come from the annotation workflow using VEP.

### `gdc-1.0.0-aliquot` Options

**Note: tumor only isn't fully supported yet**

```
usage: GDC Aliquot MAF Tools VcfToAliquotMaf gdc-1.0.0-aliquot
       [-h] [--tumor_only] [-t TUMOR_VCF_ID] [-n NORMAL_VCF_ID] --caller_id
       CALLER_ID --src_vcf_uuid SRC_VCF_UUID --case_uuid CASE_UUID
       --tumor_submitter_id TUMOR_SUBMITTER_ID --tumor_aliquot_uuid
       TUMOR_ALIQUOT_UUID --tumor_bam_uuid TUMOR_BAM_UUID
       [--normal_submitter_id NORMAL_SUBMITTER_ID]
       [--normal_aliquot_uuid NORMAL_ALIQUOT_UUID]
       [--normal_bam_uuid NORMAL_BAM_UUID] [--sequencer SEQUENCER]
       --maf_center MAF_CENTER --biotype_priority_file BIOTYPE_PRIORITY_FILE
       --effect_priority_file EFFECT_PRIORITY_FILE [--custom_enst CUSTOM_ENST]
       [--dbsnp_priority_db DBSNP_PRIORITY_DB] --reference_fasta
       REFERENCE_FASTA --reference_fasta_index REFERENCE_FASTA_INDEX
       [--reference_context_size REFERENCE_CONTEXT_SIZE]
       [--cosmic_vcf COSMIC_VCF] [--non_tcga_exac_vcf NON_TCGA_EXAC_VCF]
       [--hotspot_tsv HOTSPOT_TSV] [--exac_freq_cutoff EXAC_FREQ_CUTOFF]
       [--gdc_blacklist GDC_BLACKLIST] [--min_n_depth MIN_N_DEPTH]
       [--gdc_pon_vcf GDC_PON_VCF] [--nonexonic_intervals NONEXONIC_INTERVALS]
       [--target_intervals TARGET_INTERVALS]

optional arguments:
  -h, --help            show this help message and exit

VCF options:
  --tumor_only          Is this a tumor-only VCF?
  -t TUMOR_VCF_ID, --tumor_vcf_id TUMOR_VCF_ID
                        Name of the tumor sample in the VCF
  -n NORMAL_VCF_ID, --normal_vcf_id NORMAL_VCF_ID
                        Name of the normal sample in the VCF
  --caller_id CALLER_ID
                        Name of the caller used to detect mutations
  --src_vcf_uuid SRC_VCF_UUID
                        The UUID of the src VCF file

Sample Metadata:
  --case_uuid CASE_UUID
                        Sample case UUID
  --tumor_submitter_id TUMOR_SUBMITTER_ID
                        Tumor sample aliquot submitter ID
  --tumor_aliquot_uuid TUMOR_ALIQUOT_UUID
                        Tumor sample aliquot UUID
  --tumor_bam_uuid TUMOR_BAM_UUID
                        Tumor sample bam UUID
  --normal_submitter_id NORMAL_SUBMITTER_ID
                        Normal sample aliquot submitter ID
  --normal_aliquot_uuid NORMAL_ALIQUOT_UUID
                        Normal sample aliquot UUID
  --normal_bam_uuid NORMAL_BAM_UUID
                        Normal sample bam UUID
  --sequencer SEQUENCER
                        The sequencer used
  --maf_center MAF_CENTER
                        The sequencing center

Annotation Resources:
  --biotype_priority_file BIOTYPE_PRIORITY_FILE
                        Biotype priority JSON
  --effect_priority_file EFFECT_PRIORITY_FILE
                        Effect priority JSON
  --custom_enst CUSTOM_ENST
                        Optional custom ENST overrides
  --dbsnp_priority_db DBSNP_PRIORITY_DB
                        DBSNP priority sqlite database
  --reference_fasta REFERENCE_FASTA
                        Reference fasta file
  --reference_fasta_index REFERENCE_FASTA_INDEX
                        Reference fasta fai file
  --reference_context_size REFERENCE_CONTEXT_SIZE
                        Number of BP to add both upstream and downstream from
                        variant for reference context
  --cosmic_vcf COSMIC_VCF
                        Optional COSMIC VCF for annotating
  --non_tcga_exac_vcf NON_TCGA_EXAC_VCF
                        Optional non-TCGA ExAC VCF for annotating and
                        filtering
  --hotspot_tsv HOTSPOT_TSV
                        Optional hotspot TSV

Filtering Options:
  --exac_freq_cutoff EXAC_FREQ_CUTOFF
                        Flag variants where the allele frequency in any ExAC
                        population is great than this value as common_in_exac
                        [0.001]
  --gdc_blacklist GDC_BLACKLIST
                        The file containing the blacklist tags and tumor
                        aliquot uuids to apply them to.
  --min_n_depth MIN_N_DEPTH
                        Flag variants where normal depth is <= INT as ndp [7].
  --gdc_pon_vcf GDC_PON_VCF
                        The tabix-indexed panel of normals VCF for applying
                        the gdc pon filter
  --nonexonic_intervals NONEXONIC_INTERVALS
                        Flag variants outside of this tabix-indexed bed file
                        as NonExonic
  --target_intervals TARGET_INTERVALS
                        Flag variants outside of these tabix-indexed bed files
                        as off_target. Use one or more times.
```

### `gdc-2.0.0-aliquot` Options

**Note: we developed this to work with VEP v102**

```
usage: GDC Aliquot MAF Tools VcfToAliquotMaf gdc-2.0.0-aliquot
       [-h] [--tumor_only] [-t TUMOR_VCF_ID] [-n NORMAL_VCF_ID] --caller_id
       CALLER_ID --src_vcf_uuid SRC_VCF_UUID --case_uuid CASE_UUID
       --tumor_submitter_id TUMOR_SUBMITTER_ID --tumor_aliquot_uuid
       TUMOR_ALIQUOT_UUID --tumor_bam_uuid TUMOR_BAM_UUID
       [--normal_submitter_id NORMAL_SUBMITTER_ID]
       [--normal_aliquot_uuid NORMAL_ALIQUOT_UUID]
       [--normal_bam_uuid NORMAL_BAM_UUID] [--sequencer SEQUENCER]
       --maf_center MAF_CENTER --biotype_priority_file BIOTYPE_PRIORITY_FILE
       --effect_priority_file EFFECT_PRIORITY_FILE [--custom_enst CUSTOM_ENST]
       --reference_fasta REFERENCE_FASTA --reference_fasta_index
       REFERENCE_FASTA_INDEX [--reference_context_size REFERENCE_CONTEXT_SIZE]
       [--cosmic_vcf COSMIC_VCF] [--hotspot_tsv HOTSPOT_TSV]
       [--entrez_gene_id_json ENTREZ_GENE_ID_JSON]
       [--gnomad_noncancer_vcf GNOMAD_NONCANCER_VCF]
       [--gnomad_af_cutoff GNOMAD_AF_CUTOFF] [--gdc_blacklist GDC_BLACKLIST]
       [--min_n_depth MIN_N_DEPTH] [--gdc_pon_vcf GDC_PON_VCF]
       [--nonexonic_intervals NONEXONIC_INTERVALS]
       [--target_intervals TARGET_INTERVALS]

optional arguments:
  -h, --help            show this help message and exit

VCF options:
  --tumor_only          Is this a tumor-only VCF?
  -t TUMOR_VCF_ID, --tumor_vcf_id TUMOR_VCF_ID
                        Name of the tumor sample in the VCF
  -n NORMAL_VCF_ID, --normal_vcf_id NORMAL_VCF_ID
                        Name of the normal sample in the VCF
  --caller_id CALLER_ID
                        Name of the caller used to detect mutations
  --src_vcf_uuid SRC_VCF_UUID
                        The UUID of the src VCF file

Sample Metadata:
  --case_uuid CASE_UUID
                        Sample case UUID
  --tumor_submitter_id TUMOR_SUBMITTER_ID
                        Tumor sample aliquot submitter ID
  --tumor_aliquot_uuid TUMOR_ALIQUOT_UUID
                        Tumor sample aliquot UUID
  --tumor_bam_uuid TUMOR_BAM_UUID
                        Tumor sample bam UUID
  --normal_submitter_id NORMAL_SUBMITTER_ID
                        Normal sample aliquot submitter ID
  --normal_aliquot_uuid NORMAL_ALIQUOT_UUID
                        Normal sample aliquot UUID
  --normal_bam_uuid NORMAL_BAM_UUID
                        Normal sample bam UUID
  --sequencer SEQUENCER
                        The sequencer used
  --maf_center MAF_CENTER
                        The sequencing center

Annotation Resources:
  --biotype_priority_file BIOTYPE_PRIORITY_FILE
                        Biotype priority JSON
  --effect_priority_file EFFECT_PRIORITY_FILE
                        Effect priority JSON
  --custom_enst CUSTOM_ENST
                        Optional custom ENST overrides
  --reference_fasta REFERENCE_FASTA
                        Reference fasta file
  --reference_fasta_index REFERENCE_FASTA_INDEX
                        Reference fasta fai file
  --reference_context_size REFERENCE_CONTEXT_SIZE
                        Number of BP to add both upstream and downstream from
                        variant for reference context
  --cosmic_vcf COSMIC_VCF
                        Optional COSMIC VCF for annotating
  --hotspot_tsv HOTSPOT_TSV
                        Optional hotspot TSV
  --entrez_gene_id_json ENTREZ_GENE_ID_JSON
                        Optional map of ensembl transcript IDs and symbols to
                        entrez gene ID
  --gnomad_noncancer_vcf GNOMAD_NONCANCER_VCF
                        Path to the bgzipped and tabix-indexed non-cancer
                        gnomAD allele frequency VCF.

Filtering Options:
  --gnomad_af_cutoff GNOMAD_AF_CUTOFF
                        Flag variants where the allele frequency in any gnomAD
                        population is greater than this value as
                        common_in_gnomAD [0.001]
  --gdc_blacklist GDC_BLACKLIST
                        The file containing the blacklist tags and tumor
                        aliquot uuids to apply them to.
  --min_n_depth MIN_N_DEPTH
                        Flag variants where normal depth is <= INT as ndp [7].
  --gdc_pon_vcf GDC_PON_VCF
                        The tabix-indexed panel of normals VCF for applying
                        the gdc pon filter
  --nonexonic_intervals NONEXONIC_INTERVALS
                        Flag variants outside of this tabix-indexed bed file
                        as NonExonic
  --target_intervals TARGET_INTERVALS
                        Flag variants outside of these tabix-indexed bed files
                        as off_target. Use one or more times.
```

## Merge per-caller aliquot MAFs

Merge two or more per-caller aliquot MAFs from the same tumor/normal pair. All
implementations use the format:

```
aliquot-maf-tools MergeAliquotMafs \
    --output_maf <path to output merged MAF> \
    <schema version> <OPTIONS>
```

### `gdc-1.0.0-aliquot-merged` Options

```
usage: GDC Aliquot MAF Tools MergeAliquotMafs gdc-1.0.0-aliquot-merged
       [-h] [--mutect2 MUTECT2] [--muse MUSE] [--vardict VARDICT]
       [--varscan2 VARSCAN2] [--somaticsniper SOMATICSNIPER] [--pindel PINDEL]
       [--min_n_depth MIN_N_DEPTH]

optional arguments:
  -h, --help            show this help message and exit
  --mutect2 MUTECT2     Path to input protected MuTect2 MAF file
  --muse MUSE           Path to input protected MuSE MAF file
  --vardict VARDICT     Path to input protected VarDict MAF file
  --varscan2 VARSCAN2   Path to input protected VarScan2 MAF file
  --somaticsniper SOMATICSNIPER
                        Path to input protected SomaticSniper MAF file
  --pindel PINDEL       Path to input protected Pindel MAF file
  --min_n_depth MIN_N_DEPTH
                        Flag variants where normal depth is <= INT as ndp.
                        This is performed after averaging depths across
                        callers [7]
```

### `gdc-2.0.0-aliquot-merged` Options

**Note: currently only sets RNA columns to unknown/na, annotator not implemented yet**

```
usage: GDC Aliquot MAF Tools MergeAliquotMafs gdc-2.0.0-aliquot-merged
       [-h] [--tumor_only] [--mutect2 MUTECT2] [--muse MUSE]
       [--vardict VARDICT] [--varscan2 VARSCAN2]
       [--somaticsniper SOMATICSNIPER] [--pindel PINDEL]
       [--min_n_depth MIN_N_DEPTH]

optional arguments:
  -h, --help            show this help message and exit
  --tumor_only          If this is a tumor-only MAF
  --mutect2 MUTECT2     Path to input protected MuTect2 MAF file
  --muse MUSE           Path to input protected MuSE MAF file
  --vardict VARDICT     Path to input protected VarDict MAF file
  --varscan2 VARSCAN2   Path to input protected VarScan2 MAF file
  --somaticsniper SOMATICSNIPER
                        Path to input protected SomaticSniper MAF file
  --pindel PINDEL       Path to input protected Pindel MAF file
  --min_n_depth MIN_N_DEPTH
                        Flag variants where normal depth is <= INT as ndp.
                        This is performed after averaging depths across
                        callers [7]
```

## Masked Merged Aliquot MAF

All of the implementations use the format:

```
aliquot-maf-tools MaskMergedAliquotMaf \
    --input_maff <path to input merged MAF> \
    --output_maf <path to output filtered MAF> \
    <schema version> <OPTIONS>
```

### `gdc-1.0.0-aliquot-merged-masked` Options

```
usage: GDC Aliquot MAF Tools MaskMergedAliquotMaf gdc-1.0.0-aliquot-merged-masked
       [-h] [--tumor_only] [--reference_fasta_index REFERENCE_FASTA_INDEX]
       [--min_callers MIN_CALLERS]

optional arguments:
  -h, --help            show this help message and exit
  --tumor_only          Is this a tumor-only VCF?
  --reference_fasta_index REFERENCE_FASTA_INDEX
                        Path to the reference fasta fai file if the input MAF
                        is not sorted
  --min_callers MIN_CALLERS
                        Minimum number of callers required [2]
```

### `gdc-2.0.0-aliquot-merged-masked` Options

```
usage: GDC Aliquot MAF Tools MaskMergedAliquotMaf gdc-2.0.0-aliquot-merged-masked
       [-h] [--tumor_only] [--reference_fasta_index REFERENCE_FASTA_INDEX]
       [--min_callers MIN_CALLERS]

optional arguments:
  -h, --help            show this help message and exit
  --tumor_only          Is this a tumor-only VCF?
  --reference_fasta_index REFERENCE_FASTA_INDEX
                        Path to the reference fasta fai file if the input MAF
                        is not sorted
  --min_callers MIN_CALLERS
                        Minimum number of callers required [2]
```
