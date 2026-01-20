use clap::Parser;
use memchr::{memchr, memrchr};
use std::io::{self, BufRead, Read, Write};

#[derive(Parser)]
#[command(
    name = "fastq-i5-rc",
    version,
    about = "Rewrites FASTQ headers by reverse-complementing the i5 (Index2 / P5) barcode",
    long_about = "A fast, streaming tool to rewrite FASTQ headers by reverse-complementing the i5 (Index2 / P5) barcode, without modifying read sequences or quality scores. Headers are expected to end with the standard Illumina `:<i7>+<i5>` format."
)]
struct Args {}

/// Return the complement of a DNA base (A,C,G,T,N), preserving case.
#[inline(always)]
const fn complement_base(b: u8) -> u8 {
    // Handles A,C,G,T,N (upper/lower). Leaves other bytes unchanged.
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'n' => b'n',
        _ => b,
    }
}

/// Reverse-complement a DNA sequence in-place.
#[inline(always)]
fn reverse_complement_in_place(buf: &mut [u8]) {
    let mut i = 0;
    let mut j = buf.len();
    while i < j {
        j -= 1;
        let a = complement_base(buf[i]);
        let b = complement_base(buf[j]);
        buf[i] = b;
        buf[j] = a;
        i += 1;
    }
}

/// Reverse-complement the i5 part of a FASTQ header line in-place
/// Header is expected to start with '@' and end with ":i7+i5\n".
/// Returns an error if the header is invalid
fn rewrite_header_i5(header: &mut [u8]) -> io::Result<()> {
    if header.is_empty() || header[0] != b'@' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTQ header: does not start with '@'",
        ));
    }

    if header.last() != Some(&b'\n') {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTQ header: missing trailing newline",
        ));
    }

    // Find last ':' in the header.
    let Some(colon_index) = memrchr(b':', header) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTQ header: missing ':' before index field",
        ));
    };
    let after_colon_index = colon_index + 1;

    // Find '+' after that last ':'.
    let Some(relative_plus_index) = memchr(b'+', &header[after_colon_index..]) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid FASTQ header: missing '+' in index field",
        ));
    };
    let after_plus_index = after_colon_index + relative_plus_index + 1;

    // i5 header is everything after '+' excluding the final newline character
    let stop_h5_index = header.len() - 1; // exclude final '\n'
    let i5 = &mut header[after_plus_index..stop_h5_index];
    reverse_complement_in_place(i5);
    Ok(())
}

/// Read one line (including the trailing '\n' if present) into line, erasing previous contents.
/// Returns number of bytes read (0 on EOF).
#[inline(always)]
fn read_line<R: Read>(reader: &mut io::BufReader<R>, line: &mut Vec<u8>) -> io::Result<usize> {
    line.clear();
    let mut total = 0usize;
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(total);
        }
        if let Some(pos) = memchr(b'\n', available) {
            // include newline
            line.extend_from_slice(&available[..=pos]);
            let consume = pos + 1;
            reader.consume(consume);
            total += consume;
            return Ok(total);
        } else {
            // consume all
            line.extend_from_slice(available);
            let consume = available.len();
            reader.consume(consume);
            total += consume;
        }
    }
}

fn main() -> io::Result<()> {
    let _args = Args::parse();
    const IO_BUFFER_BYTES: usize = 64 * 1024; // 64 kB buffer for I/O
    let stdin = io::stdin();
    let mut input = io::BufReader::with_capacity(IO_BUFFER_BYTES, stdin.lock());
    let stdout = io::stdout();
    let mut output = io::BufWriter::with_capacity(IO_BUFFER_BYTES, stdout.lock());

    // Buffer for a FASTQ record line (a record is 4 lines where the first line is the header)
    let mut line = Vec::<u8>::with_capacity(1);
    const N_LINES_PER_RECORD: usize = 4;

    loop {
        // rewrite header line
        if read_line(&mut input, &mut line)? == 0 {
            break; // no header: EOF
        }
        rewrite_header_i5(&mut line)?;
        output.write_all(&line)?;

        // copy remaining lines of the FASTQ record unchanged
        for _ in 1..N_LINES_PER_RECORD {
            if read_line(&mut input, &mut line)? == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated FASTQ record (expected 4 lines)",
                ));
            }
            output.write_all(&line)?;
        }
    }

    output.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::actacttgag(
        b"@VH00821:6:AACCCKLM5:1:1101:18231:1000 1:N:0:TCTTGAGGTT+ACTACTTGAG\n",
        b"@VH00821:6:AACCCKLM5:1:1101:18231:1000 1:N:0:TCTTGAGGTT+CTCAAGTAGT\n"
    )]
    #[case::a(b"@r1 1:N:0:AAAA+A\n", b"@r1 1:N:0:AAAA+T\n")]
    #[case::ac(b"@r1 1:N:0:AAAA+AC\n", b"@r1 1:N:0:AAAA+GT\n")]
    #[case::acg(b"@r1 1:N:0:AAAA+ACG\n", b"@r1 1:N:0:AAAA+CGT\n")]
    #[case::acgt(b"@r1 1:N:0:AAAA+ACGT\n", b"@r1 1:N:0:AAAA+ACGT\n")]
    #[case::acgt_lowercase(b"@r2 1:N:0:CCCC+acgt\n", b"@r2 1:N:0:CCCC+acgt\n")]
    #[case::nnnn(b"@r3 1:N:0:GGGG+NNNN\n", b"@r3 1:N:0:GGGG+NNNN\n")]
    #[case::actg_mixedcase(b"@r4 1:N:0:TTTT+AcTg\n", b"@r4 1:N:0:TTTT+cAgT\n")]
    #[case::extra_colons(
        b"@inst:run:flow:lane:tile:x:y 1:N:0:AAAA+TTTT\n",
        b"@inst:run:flow:lane:tile:x:y 1:N:0:AAAA+AAAA\n"
    )]
    #[case::empty_i5(b"@pyt1 1:N:0:AAAA+\n", b"@pyt1 1:N:0:AAAA+\n")]
    #[case::single_a(b"@pyt2 1:N:0:AAAA+A\n", b"@pyt2 1:N:0:AAAA+T\n")]
    #[case::single_n(b"@pyt3 1:N:0:AAAA+N\n", b"@pyt3 1:N:0:AAAA+N\n")]
    #[case::mixed_case_short(b"@pyt4 1:N:0:AAAA+AaCg\n", b"@pyt4 1:N:0:AAAA+cGtT\n")]
    #[case::acgtn(b"@pyt5 1:N:0:AAAA+AcgTN\n", b"@pyt5 1:N:0:AAAA+NAcgT\n")]
    #[case::all_as(b"@pyt6 1:N:0:AAAA+AAAA\n", b"@pyt6 1:N:0:AAAA+TTTT\n")]
    #[case::all_cs(b"@pyt7 1:N:0:AAAA+CCCC\n", b"@pyt7 1:N:0:AAAA+GGGG\n")]
    #[case::at_repeat(b"@pyt8 1:N:0:AAAA+ATATAT\n", b"@pyt8 1:N:0:AAAA+ATATAT\n")]
    #[case::cg_repeat(b"@pyt9 1:N:0:AAAA+CGCGCG\n", b"@pyt9 1:N:0:AAAA+CGCGCG\n")]
    #[case::ns_flanking(b"@pyt10 1:N:0:AAAA+NNACGTNN\n", b"@pyt10 1:N:0:AAAA+NNACGTNN\n")]
    #[case::general_atcacg(b"@pyt11 1:N:0:AAAA+ATCACG\n", b"@pyt11 1:N:0:AAAA+CGTGAT\n")]
    #[case::general_ttaggc(b"@pyt12 1:N:0:AAAA+TTAGGC\n", b"@pyt12 1:N:0:AAAA+GCCTAA\n")]
    fn rewrite_header_i5_valid(
        #[case] input: &[u8],
        #[case] expected: &[u8],
    ) -> std::io::Result<()> {
        let mut header = input.to_vec();
        rewrite_header_i5(&mut header)?;
        assert_eq!(
            String::from_utf8_lossy(&header),
            String::from_utf8_lossy(expected),
            "\n\n     test input: {}\n    test output: {}\nexpected output: {}\n",
            String::from_utf8_lossy(input),
            String::from_utf8_lossy(&header),
            String::from_utf8_lossy(expected),
        );
        Ok(())
    }

    #[rstest]
    #[case::empty_header(b"@\n", "missing ':'")]
    #[case::no_colon(b"@r6 no_index_here\n", "missing ':'")]
    #[case::no_plus(b"@r5 1:N:0:AAAA\n", "missing '+'")]
    #[case::no_newline(b"@r7 1:N:0:CCCC+AGTC", "missing trailing newline")]
    fn rewrite_header_i5_invalid(#[case] input: &[u8], #[case] msg_substr: &str) {
        let mut header = input.to_vec();
        let err = rewrite_header_i5(&mut header).expect_err("expected rewrite_header_i5 to fail");
        let msg = err.to_string();
        assert!(
            msg.contains(msg_substr),
            "error message did not contain expected substring.\n  expected: {msg_substr}\n  got: {msg}"
        );
    }
}
