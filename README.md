# nf-run-snpAD

## Overview
`nf-run-snpAD` is a Nextflow pipeline designed to run SNP analysis on BAM files. This pipeline processes samples with or without UDG treatment based on the provided input.

## Prerequisites
Nextflow: Ensure you have Nextflow installed. You can install it by running:

```bash
curl -s https://get.nextflow.io | bash
```
- Java: Nextflow requires Java. Make sure you have Java (version 8 or higher) installed on your system.
- Modules: This pipeline includes several local modules that must be available in the designated paths.

## Installation
Clone the repository or download the pipeline files:

```bash
git clone https://github.com/jbv2/run-snpAD.git
cd nf-run-snpAD
```

## Usage
To run the pipeline, use the following command:

```bash
nextflow run main.nf -profile <profile_name> --input_tsv <path_to_input_tsv>
```

### Parameters
`--input_tsv`: Path to the TSV file containing sample information (sample_id, inputbam, bai, udg).
`-profile`: Specify a Nextflow profile (e.g., test, docker, etc.) for execution.

> Input TSV Format
> Your input TSV file should have the following format:
> ```bash
> sample_id	inputbam	bai	udg
> sample1	sample1.bam	sample1.bai	true
> sample2	sample2.bam	sample2.bai	false
> ```

## Example Command
```bash
nextflow run main.nf -profile test --input_tsv test/data/test.tsv
```

## Output
The pipeline generates the following outputs:

* VCF files containing the called variants for each sample, stored in the defined output directory.

### Note
Ensure that your input files and paths are correctly set, and that you have the necessary permissions to read the files and write the outputs.

## Support
For issues or questions, please open an issue on the GitHub repository.

## License
This project is licensed under the MIT License - see the LICENSE file for details.