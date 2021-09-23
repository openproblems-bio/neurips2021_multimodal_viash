# $par_task_name - Starter Kit for $par_language_name Users

Full documentation for the competition, including much of the information here, can be found online 
at [openproblems.bio/neurips_docs/](https://openproblems.bio/neurips_docs/). The documentation for 
Viash is available at [viash.io/docs](https://viash.io/docs).
​
## Getting Started
​
Check the [Quickstart](https://openproblems.bio/neurips_docs/submission/quickstart/) to create and upload your first submission to EvalAI.
​
Check the following links for more information:
​
- [Starter kit contents](https://openproblems.bio/neurips_docs/submission/starter_kit_contents/)
- [Development process](https://openproblems.bio/neurips_docs/submission/development_process/)
- [Submit to EvalAI](https://eval.ai/web/challenges/challenge-page/1111/submission)
​
## Folder Structure
​
```
├── LICENSE                                 # MIT License
├── README.md                               # Some starter information
├── bin/                                    # Binaries needed to generate a submission
│   ├── check_format
│   ├── nextflow
│   └── viash
├── config.vsh.yaml                         # Viash configuration file
├── script.py                               # Script containing your method
├── sample_data/                            # Small sample datasets for unit testing and debugging
│   ├── openproblems_bmmc_cite_starter/     # Contains H5AD files for CITE data
│   └── openproblems_bmmc_multiome_starter/ # Contains H5AD files for multiome data
├── scripts/                                # Scripts to test, generate, and evaluate a submission
│   ├── 0_sys_checks.sh                     # Checks that necessary software installed
│   ├── 1_unit_test.sh                      # Runs the unit tests in test.py
│   ├── 2_generate_submission.sh            # Generates a submission pkg by running your method on validation data
│   ├── 3_evaluate_submission.sh            # (Optional) Scores your method locally
│   └── nextflow.config                     # Configurations for running Nextflow locally
└── test.py                                 # Default unit tests. Feel free to add more tests, but don't remove any.
```

