# ExpansionHunterDenovo LRDN

A tool for local read depth normalization of ExpansionHunter denovo repeat expansion calls.

## Requirements

There are several Python packages required to run the `ehdn_lrdn.py` and `merge_lrdn_output.py` scripts. These can easily be installed using `pip` and the `requirements.txt` file provided:

```
$ pip install -r requirements.txt
```

## Usage

First, run the `ehdn_lrdn.py` script on each sample:

```
$ python ehdn_lrdn.py case_sample1.bam case_sample1 ehdn_output.tsv case_sample1.str_profile.json lrdn_output_dir/
```

Then, run the merge script:

```
$ python merge_lrdn_output.py lrdn_output_dir/ lrdn_output.tsv
```
