use assert_cmd::cargo::*;

#[test]
fn valid() {
    let input = b"@r1 1:N:0:AAAA+ACTACTTGAG\n\
ACGT\n\
+\n\
!!!!\n\
@r2 1:N:0:CCCC+atcacg\n\
TGCA\n\
+\n\
####\n";

    let expected = b"@r1 1:N:0:AAAA+CTCAAGTAGT\n\
ACGT\n\
+\n\
!!!!\n\
@r2 1:N:0:CCCC+cgtgat\n\
TGCA\n\
+\n\
####\n";

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");
    cmd.write_stdin(input)
        .assert()
        .success()
        .stdout(&expected[..]);
}

#[test]
fn invalid() {
    // Missing '+' after the last ':' in the header
    let input = b"@r1 1:N:0:AAAAACGT\n\
ACGT\n\
+\n\
!!!!\n";

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");
    cmd.write_stdin(input)
        .assert()
        .failure()
        .stderr(predicates::str::contains("missing '+'"));
}
