# LUAD

Raw data in `raw` have been downloaded from the following links: 

- `expression.csv` from `https://grand.networkmedicine.org/tissues/Colon_transverse_tissue/` [accessed April 6 2023]
- `motif.txt` from `https://grand.networkmedicine.org/tissues/Colon_transverse_tissue/` [accessed April 6 2023] 
- `TF-Target-information.txt` from `http://bioinfo.life.hust.edu.cn/hTFtarget/#!/` [accessed April 5 2023]


TF-gene priors are from CIS-BP, and PPI has been extracted from STRING. 

The other files are computed either by `preprocessing.py`, by the corresponding R file, or by the jupyter notebook using the option `cache=False`. 