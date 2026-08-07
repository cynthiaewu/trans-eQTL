# Trex-QTL

This repository contains the following utilities for simulating expression datasets and running various methods for trans-eQTL detection:

* [`trexqtl-simulate`](simulate/README.md): Tool for simulating expression datasets for a cohort of individuals
* [`trexqtl-run`](trexqtl/README.md): Tool for detecting trans-eQTL using either the published CPMA method or our new mixture model based method (Trex-QTL). 

# Install

To install Trex-QTL locally, run:

```
python3 setup.py install --user
```

This will install the commands `trexqtl-simulate` and `trexqtl-run`. Type `trexqtl-simulate --help` or `trexqtl-run --help` for usage information.
