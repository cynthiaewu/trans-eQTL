# trans-eQTL

** Under construction **

This repository contains the following utilities for simulating expression datasets and running various methods for trans-eQTL detection:

* `xQTL-simulate`: Tool for simulating expression datasets for a cohort of individuals
* `xQTL-run`: Tool for detecting trans-eQTL using either the published CPMA method or our new mixture model based method (xQTL). (**Not yet implemented**)

# Install

To install xQTL locally, run:

```
python3 setup.py install --user
```

This will install the commands `xQTL-simulate` and `xQTL-run`.

# Developer notes

To run tests:

```
python -m pytest --cov=. --cov-report term-missing
```
