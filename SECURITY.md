# Security policy

Before the first tagged release, security fixes target the `main` branch.
After publication, the latest tagged release is supported.

Use **Report a vulnerability** on the LOCUSKIT repository's Security tab if
private reporting is available. Otherwise, contact a LOCUSKIT maintainer
privately. Include the package version and commit, Python and tool versions,
a minimal reproduction, and the expected impact when possible. Keep sensitive
genomic data and credentials out of public issues and logs.

The base LOCUSKIT installation has no third-party runtime dependencies.
The optional `compare` extra installs JCVI and its scientific dependency tree;
use `uv sync --locked --extra all` for the reviewed dependency set and keep the
lock file updated. As reviewed on 2026-08-31, JCVI 1.6.6 pulls `diskcache 5.6.3`
through an unused `genomepy` dependency.
[CVE-2025-69872 / GHSA-w8v5-vhqr-4h9v](https://github.com/advisories/GHSA-w8v5-vhqr-4h9v)
has no fixed release and requires
an attacker to have write access to the cache directory; kfFractBias does not
import or create a DiskCache. CI allows only that advisory explicitly and
fails for any additional known vulnerability. Its pip-audit identifier is
`PYSEC-2026-2447`. Recheck the advisory when refreshing dependencies; the review
date is not a promise that upstream status will remain unchanged.

External FASTA, GFF, BED, and synteny files are untrusted inputs and should
not share a writable output path.
