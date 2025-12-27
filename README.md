# NCAA Token v1.1 Pipeline (Lite)

Directories:
- input/    : immutable raw inputs + manifest
- config/   : token v1.1 rules, thresholds, vocab, fixed pH settings
- shards/   : per-SID-range outputs (params tarball + token parquet + validate report)
- registry/ : merged final registry (v1.1)
- qc/       : token distribution & sanity checks
- logs/     : per-run logs
