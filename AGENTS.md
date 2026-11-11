<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=12; sha256=d38271eaba87f3f2c5c92b6c4ea821eca662d38f6f0db388fda5270edb5e9afa -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Prefer the simplest design that meets established requirements. Avoid
  complexity for speculative needs; extract helpers when they improve clarity
  or enable reuse.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Validate untrusted input at the responsible boundary; rely on established
  invariants internally. Revalidate when mutation or a new trust boundary
  invalidates those guarantees.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- Optimize bottlenecks in RAM usage or CPU time whenever found. For software
  performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->

# Working in LOCUSKIT

## Start here

Read [CONTRIBUTING.md](CONTRIBUTING.md) for the development environment and
check policy, then the relevant section of [README.md](README.md) for the CLI.
For numerical changes read [methods](docs/methods.md); for parsing or output
changes read [formats](docs/formats.md) and [migration](docs/migration.md).
Use [the documentation index](docs/README.md) for other topics, not the upstream
CoGe web-tool instructions.

The LOCUSKIT CLI starts at `src/locuskit/cli.py`. Read [the neighborhood format](docs/neighborhoods.md) for locus or context changes.
The FractBias execution path starts at `src/kffractbias/cli.py`; `analysis.py` orchestrates
calculation and `profiles.py` produces retention rows. `io.py` owns input
validation, `run.py` owns locking/staging/rollback, and `jcvi.py` delegates
self-synteny to `selfscan.py`. Start with the module relevant to the change.

## Run and verify

Work from the repository root, with Python 3.11+ on Linux/macOS and a local
POSIX filesystem. Initial setup: `uv sync --locked --extra test --extra plot`.
Keep an existing comparison environment with `--extra all` instead of `plot`;
syncing fewer extras removes optional packages. Subsequent checks use
`uv run --no-sync` to preserve the prepared environment.

- Standard check before push: `uv run --no-sync python scripts/check.py`.
- Small CLI execution with temporary outputs:
  `uv run --no-sync python -m pytest -q -rs tests/test_documentation.py`.
- Lint: `uv run --no-sync ruff check src tests scripts`.
- Format check: `uv run --no-sync ruff format --check src tests scripts`.
- Types: `uv run --no-sync mypy src`.

Use the [change-to-test table](CONTRIBUTING.md#choosing-checks) to select focused
tests and decide whether distribution or real-alignment checks are needed.
The standard check already runs lint, format, types, and the fast test suite;
do not repeat them after it passes without a reason. Skipped tests do not
verify optional functionality. For CLI/documentation changes, the local
[check-cli-docs skill](.agents/skills/check-cli-docs/SKILL.md) describes the
existing example checks and their limits.

## Preserve scientific and compatibility contracts

- Preserve gene-based complete windows, denominator filtering before ranking,
  binary retention per query sequence, and deterministic dense TSV rows.
  BED coordinates and output ranks use different conventions; see methods.
- Quota is target:query. Self retention has a different interpretation from
  outgroup fractionation; its diagonal filter is chromosome-local and its
  depth constraint combines both arms. Do not tune biological defaults to
  make a test pass or infer a user's quota/depth.
- Preserve CLI options, TSV columns/order/precision, summary schema and hashes,
  and the Python `AnalysisConfig`/`AnalysisResult` behavior (library collection
  defaults differ from CLI streaming). Intentional compatibility changes need
  corresponding tests and updates to formats/migration documentation.
- Preserve input snapshots, per-prefix locking, and rollback. Never bypass
  an active lock; file existence alone does not indicate lock ownership.

## Boundaries and completion

Keep existing user edits. Do not edit biological inputs, retained synteny work,
result/staging/lock files, local environments, caches, `dist/`, or `build/` as
source. Run examples through the temporary-directory tests; standalone CLI
examples can replace results with the same prefix. Do not introduce real data
downloads or long analyses into routine verification. Change dependency locks,
versions, or research parameters only when the task calls for it, following
CONTRIBUTING for version changes.

Before finishing, inspect the diff for unintended files and run the required
checks. Report changed behavior/files, commands and outcomes, skipped checks
and missing prerequisites, and any compatibility implications. Distinguish
executed checks from static inspection; do not claim local results cover the
entire CI platform matrix.
