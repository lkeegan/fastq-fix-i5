# fastq-fix-i5

A fast, streaming tool to **rewrite FASTQ headers by
reverse-complementing the i5 (Index2 / P5) barcode**,
without modifying read sequences or quality scores.

Headers are expected to end with the standard Illumina
`:<i7>+<i5>` format.

This tool can be useful when mixing FASTQs from different sequencing
platforms (e.g. Illumina and AVITI) where i5 orientation conventions differ.

## What this tool does

- Reads FASTQ from **stdin**
- Writes FASTQ to **stdout**
- For each record:
    - parses the FASTQ header
    - finds the final `:<i7>+<i5>` field
    - **reverse-complements only the i5 part**
- Leaves everything else unchanged

## Installation

Download a pre-compiled binary for your platform.

## Use

To run `fastq-fix-i5`, simply pipe your FASTQ data into it:

```bash
fastq-fix-i5 < input.fastq > output.fastq
```

If your input FASTQ is compressed, on linux you can use `pigz` to
stream the data through `fastq-fix-i5`:

```bash
pigz -dc input.fastq.gz | fastq-fix-i5 | pigz -c > output.fastq.gz
```

## Performance

`fastq-fix-i5` is designed to be fast and memory-efficient.
It processes FASTQ data in a streaming fashion,
which allows it to process millions of reads per second
using a single CPU core and <5MB of memory.
