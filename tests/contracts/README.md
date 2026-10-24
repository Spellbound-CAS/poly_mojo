# API declaration fixtures

These JSON snapshots are inputs to `tests/test_contract_p.py`. They record
current exact and numerical declarations and historical baselines used to
detect unintended source-compatibility changes. They are test data and are
not published with the documentation.

Run from the repository root:

```sh
pixi run --locked python tests/test_contract_p.py
```

The test has no snapshot-update mode. Change a snapshot only with a reviewed
API change and an accompanying consumer migration. See
[contributing](../../CONTRIBUTING.md) and the
[public compatibility guide](../../docs/content/reference/contract-p.md).
