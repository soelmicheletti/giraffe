# LUAD

Raw data in `raw` have been downloaded from the following paths in [Netbooks](http://netbooks.networkmedicine.org/hub/login) (see the also the [GitHub repo](https://github.com/netZoo/netbooks)):

- `expression.csv` from `/opt/data/netZooPy/krcc/scaled_gene_expression_KRCC.csv` [accessed April 5 2023]
- `motif.txt` from `https://grand.networkmedicine.org/tissues/Kidney_cortex_tissue/motif.txt` [accessed April 5 2023] 


TF-gene priors are from CIS-BP, and PPI has been extracted from STRING. 

The files `expr.csv`, `expr_t.csv` and `prior.csv` are copies of expression and prior reformatted for convenience of R scripts.
`R_tigress.csv` is computed via `tigress.r`. The other files are computed either by `preprocessing.py` or by the jupyter notebook using the option `cache=False`. 