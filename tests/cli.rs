use assert_cmd::cargo::*;

#[test]
fn valid_two_records() {
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
    // piping the output in again should recover the original input
    cmd.write_stdin(expected)
        .assert()
        .success()
        .stdout(&input[..]);
}

#[test]
fn valid_5k_records_roundtrip() {
    let record = b"@r1 1:N:0:AAAA+ACTACTTGAG\n\
ACGT\n\
+\n\
!!!!\n";

    let input = record.repeat(5000);

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");

    // First pass: transform the input
    let output = cmd
        .write_stdin(input.clone())
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    // Second pass: piping the output back in should recover the original input
    cmd.write_stdin(output).assert().success().stdout(input);
}

#[test]
fn invalid_missing_plus() {
    // Missing '+' after the last ':' in the header
    let input = b"@r1 1:N:0:AAAAACGT\n\
ACGT\n\
+\n\
!!!!\n";

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");
    cmd.write_stdin(input)
        .assert()
        .failure()
        .stderr(predicates::str::contains("'+'"));
}

#[test]
fn invalid_truncated_record() {
    let input = b"@r1 1:N:0:CCCC+atcacg\n\
ACGT\n\
...\n";

    let mut cmd = cargo_bin_cmd!("fastq-fix-i5");
    cmd.write_stdin(input)
        .assert()
        .failure()
        .stderr(predicates::str::contains("truncated"));
}
