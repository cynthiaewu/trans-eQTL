# xQTL


This repository contains the following utilities for simulating expression datasets and running various methods for trans-eQTL detection:

* [`xQTL-simulate`](simulate/README.md): Tool for simulating expression datasets for a cohort of individuals
* [`xQTL-run`](xQTL/README.md): Tool for detecting trans-eQTL using either the published CPMA method or our new mixture model based method (xQTL). 

# Install

To install xQTL locally, run:

```
python3 setup.py install --user
```

This will install the commands `xQTL-simulate` and `xQTL-run`. Type `xQTL-simulate --help` or `xQTL-run --help` for usage information.

# Developer notes

To run tests:

```
python -m pytest --cov=. --cov-report term-missing
```
