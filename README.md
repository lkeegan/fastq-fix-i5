# fastq-fix-i5

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/fastq-fix-i5/README.html)
[![Crates.io Version](https://img.shields.io/crates/v/fastq-fix-i5)](https://crates.io/crates/fastq-fix-i5)
[![codecov](https://codecov.io/gh/ssciwr/fastq-fix-i5/graph/badge.svg?token=xn3GEKWMO1)](https://codecov.io/gh/ssciwr/fastq-fix-i5)

A fast, streaming tool to **rewrite FASTQ headers by
reverse-complementing the i5 (Index2 / P5) barcode**,
without modifying read sequences or quality scores.

Headers are expected to end with the standard Illumina
`:<i7>+<i5>` format.

This tool can be useful when mixing FASTQs from different sequencing
platforms (e.g. Illumina and AVITI) where i5 orientation conventions differ.

## What this tool does

- Reads FASTQ from **stdin**
- For each record:
    - parses the FASTQ header
    - finds the final `:<i7>+<i5>` field
    - **reverse-complements only the i5 part**
    - leaves everything else unchanged
- Writes FASTQ to **stdout**

If it encounters a header that does not conform to the expected format,
or a truncated FASTQ record,
it will print an error message to stderr and exit with a non-zero status code.

## Installation

To install from bioconda:

```
conda install -c bioconda fastq-fix-i5
```

To build and install using the rust package manager cargo:

```
cargo install fastq-fix-i5
```

Alternatively you can download a [pre-compiled binary](https://github.com/ssciwr/fastq-fix-i5/releases/latest) for your platform.

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
which allows it to process millions of reads per second using a single CPU core,
using a tiny fixed amount (< 3MB) of RAM regardless of input size.

Here are some example benchmarks for processing 10 million synthetic fastq reads on a standard desktop computer:

| Command                                                           | Time | Throughput           |
|-------------------------------------------------------------------|------|----------------------|
| `fastq-fix-i5 < 10m.fastq > /dev/null`                            | 3.2s | 3.1 million reads/s  |
| `fastq-fix-i5 < 10m.fastq > out.fastq`                            | 9.2s | 1.1 million reads/s  |
| `pigz -dc 10m.fastq.gz \| fastq-fix-i5 > /dev/null`               | 4.0s | 2.5 million reads/s  |
| `pigz -dc 10m.fastq.gz \| fastq-fix-i5 \| pigz -c > out.fastq.gz` | 6.5s | 1.5 million reads/s  |

In addition, since the runtime in these benchmarks is largely I/O-bound,
inserting `fastq-fix-i5` into an existing Unix pipe (`|`) between other commands
will have minimal impact on the overall speed:

```bash
Without fastq-fix-i5: pigz -dc 10m.fastq.gz | pigz -c > out.fastq.gz
  Time (mean ± σ):      6.445 s ±  0.079 s    [User: 39.832 s, System: 17.983 s]
  Range (min … max):    6.327 s …  6.582 s    10 runs

With fastq-fix-i5: pigz -dc 10m.fastq.gz | fastq-fix-i5 | pigz -c > out.fastq.gz
  Time (mean ± σ):      6.526 s ±  0.014 s    [User: 49.848 s, System: 16.374 s]
  Range (min … max):    6.508 s …  6.543 s    10 runs
````

The synthetic FASTQ file used above was generated using the following command:

```bash
yes "@VH00821:6:AACCCKLM5:1:1101:18231:1000 1:N:0:TCTTGAGGTT+ATCACGATCA\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" | head -n $((10000000 * 4)) > 10m.fastq
```
