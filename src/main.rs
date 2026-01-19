use std::io::{self, Read, Write, BufRead};
use memchr::{memchr,memrchr};

/// Return the complement of a DNA base (A,C,G,T,N), preserving case.
#[inline(always)]
fn complement_base(b: u8) -> u8 {
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

/// Rewrite a FASTQ header line expected to end with "...:i7+i5\n" (or without final newline).
/// We reverse-complement i5 only and write the modified header to `out`.
/// Returns false if pattern not found (no rewrite); true if rewritten.
fn rewrite_header_i5(header: &mut Vec<u8>) -> bool {
    if header.is_empty() || header[0] != b'@' {
        // Not a FASTQ header; pass through unchanged.
        return false;
    }

    // Check for trailing newline.
    let has_nl = header.last() == Some(&b'\n');

    // Find last ':' in the header.
    let Some(j) = memrchr(b':', header) else {
        return false;
    };

    // Find '+' after that last ':'.
    let Some(rel_k) = memchr(b'+', &header[j + 1..]) else {
        return false;
    };

    // i5 header is everything after '+' excluding the newline if present
    let i5_start = j + 1 + rel_k + 1;
    let i5_end = if has_nl {
        header.len() - 1
    } else {
        header.len()
    };
    let i5 = &mut header[i5_start..i5_end];
    reverse_complement_in_place(i5);
    true
}

fn main() -> io::Result<()> {
    let read_buf_size = 64*1024;  // 1 MiB buffer for I/O
    let write_buf_size = 64*1024;  // 1 MiB buffer for I/O
    let stdin = io::stdin();
    let mut input = io::BufReader::with_capacity(read_buf_size, stdin.lock());
    let stdout = io::stdout();
    let mut output = io::BufWriter::with_capacity(write_buf_size, stdout.lock());

    // Buffers for the 4 FASTQ record lines that constitute one read
    let mut h = Vec::<u8>::with_capacity(256);
    let mut s = Vec::<u8>::with_capacity(256);
    let mut p = Vec::<u8>::with_capacity(256);
    let mut q = Vec::<u8>::with_capacity(256);

    loop {
        h.clear();
        let n = read_line(&mut input, &mut h)?;
        if n == 0 {
            break; // EOF
        }

        s.clear();
        p.clear();
        q.clear();

        if read_line(&mut input, &mut s)? == 0
            || read_line(&mut input, &mut p)? == 0
            || read_line(&mut input, &mut q)? == 0
        {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ record (expected 4 lines)",
            ));
        }

        rewrite_header_i5(&mut h);
        output.write_all(&h)?;
        output.write_all(&s)?;
        output.write_all(&p)?;
        output.write_all(&q)?;
    }

    output.flush()?;
    Ok(())
}

/// Read one line (including the trailing '\n' if present) into buf.
/// Returns number of bytes read (0 on EOF).
#[inline(always)]
fn read_line<R: Read>(r: &mut io::BufReader<R>, buf: &mut Vec<u8>) -> io::Result<usize> {
    buf.clear();
    let mut total = 0usize;
    loop {
        let available = r.fill_buf()?;
        if available.is_empty() {
            return Ok(total);
        }
        if let Some(pos) = memchr(b'\n', available) {
            // include newline
            buf.extend_from_slice(&available[..=pos]);
            let consume = pos + 1;
            r.consume(consume);
            total += consume;
            return Ok(total);
        } else {
            // consume all
            buf.extend_from_slice(available);
            let consume = available.len();
            r.consume(consume);
            total += consume;
        }
    }
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
    #[case::acgt(
        b"@r1 1:N:0:AAAA+ACGT\n",
        b"@r1 1:N:0:AAAA+ACGT\n"
    )]
    #[case::acgt_lowercase(
        b"@r2 1:N:0:CCCC+acgt\n",
        b"@r2 1:N:0:CCCC+acgt\n"
    )]
    #[case::nnnn(
        b"@r3 1:N:0:GGGG+NNNN\n",
        b"@r3 1:N:0:GGGG+NNNN\n"
    )]
    #[case::actg_mixedcase(
        b"@r4 1:N:0:TTTT+AcTg\n",
        b"@r4 1:N:0:TTTT+cAgT\n"
    )]
    #[case::extra_colons(
        b"@inst:run:flow:lane:tile:x:y 1:N:0:AAAA+TTTT\n",
        b"@inst:run:flow:lane:tile:x:y 1:N:0:AAAA+AAAA\n"
    )]
    #[case::no_plus(
        b"@r5 1:N:0:AAAA\n",
        b"@r5 1:N:0:AAAA\n"
    )]
    #[case::no_colon(
        b"@r6 no_index_here\n",
        b"@r6 no_index_here\n"
    )]
    #[case::no_newline(
        b"@r7 1:N:0:CCCC+AGTC",
        b"@r7 1:N:0:CCCC+GACT"
    )]
    #[case::empty_header(
        b"@\n",
        b"@\n"
    )]
    #[case::empty_i5(
        b"@pyt1 1:N:0:AAAA+\n",
        b"@pyt1 1:N:0:AAAA+\n"
    )]
    #[case::single_a(
        b"@pyt2 1:N:0:AAAA+A\n",
        b"@pyt2 1:N:0:AAAA+T\n"
    )]
    #[case::single_n(
        b"@pyt3 1:N:0:AAAA+N\n",
        b"@pyt3 1:N:0:AAAA+N\n"
    )]
    #[case::mixed_case_short(
        b"@pyt4 1:N:0:AAAA+AaCg\n",
        b"@pyt4 1:N:0:AAAA+cGtT\n"
    )]
    #[case::acgtn(
        b"@pyt5 1:N:0:AAAA+AcgTN\n",
        b"@pyt5 1:N:0:AAAA+NAcgT\n"
    )]
    #[case::all_as(
        b"@pyt6 1:N:0:AAAA+AAAA\n",
        b"@pyt6 1:N:0:AAAA+TTTT\n"
    )]
    #[case::all_cs(
        b"@pyt7 1:N:0:AAAA+CCCC\n",
        b"@pyt7 1:N:0:AAAA+GGGG\n"
    )]
    #[case::at_repeat(
        b"@pyt8 1:N:0:AAAA+ATATAT\n",
        b"@pyt8 1:N:0:AAAA+ATATAT\n"
    )]
    #[case::cg_repeat(
        b"@pyt9 1:N:0:AAAA+CGCGCG\n",
        b"@pyt9 1:N:0:AAAA+CGCGCG\n"
    )]
    #[case::ns_flanking(
        b"@pyt10 1:N:0:AAAA+NNACGTNN\n",
        b"@pyt10 1:N:0:AAAA+NNACGTNN\n"
    )]
    #[case::general_atcacg(
        b"@pyt11 1:N:0:AAAA+ATCACG\n",
        b"@pyt11 1:N:0:AAAA+CGTGAT\n"
    )]
    #[case::general_ttaggc(
        b"@pyt12 1:N:0:AAAA+TTAGGC\n",
        b"@pyt12 1:N:0:AAAA+GCCTAA\n"
    )]
    fn rewrite_header_i5_cases(#[case] input: &[u8], #[case] expected: &[u8]) {
        let mut header = input.to_vec();
        rewrite_header_i5(&mut header);
        assert_eq!(
            header.as_slice(),
            expected,
            "input = {:?}",
            std::str::from_utf8(input).unwrap_or("<non-utf8>")
        );
    }
}
