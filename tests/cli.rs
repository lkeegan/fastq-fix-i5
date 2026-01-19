use assert_cmd::cargo::*;

#[test]
fn two_valid_reads() {
    let input = b"@r1 1:N:0:AAAA+ACGT\n\
ACGT\n\
+\n\
!!!!\n\
@r2 1:N:0:CCCC+nnnn\n\
TGCA\n\
+\n\
####\n";

    // i5 RC:
    // ACGT -> ACGT (RC of ACGT is ACGT)
    // nnnn -> nnnn
    let expected = b"@r1 1:N:0:AAAA+ACGT\n\
ACGT\n\
+\n\
!!!!\n\
@r2 1:N:0:CCCC+nnnn\n\
TGCA\n\
+\n\
####\n";

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");
    cmd.write_stdin(input)
        .assert()
        .success()
        .stdout(&expected[..]);
}
