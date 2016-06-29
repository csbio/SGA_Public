Key for standard SGA file formats. 
===================================

Newer formats should have a header. Several alternative
formats were 'standard' for a short while, but didn't 
survive and are not described here.

Conveniently, these all have a unique number of columns
* 9  Cols: Database dump, input to compute_sgascore.m
* 11 Cols: Final format of 2016 paper
* 12 Cols: Output of compute_sgascore.m
* 13 Cols: Final format of 2010 paper


INPUT to compute_sgascore.m
---------------------------
(rawdata from database)

### 9 col 
1. Query Orf
2. Array Orf
3. Array plate num (e.g. 354)
4. Set ID
5. Unique plate ID
6. Batch ID
7. row    (1 : 32)
8. column (1 : 48)
9. colony_size (in pixels)

OUTPUT of compute_sgascore.m 
-----------------------------
(scored raw data)

### 12 col
1. Query Orf
2. Array Orf
3. e_score\*
4. e_score (std)
5. p-value
6. Query smf
7. Query smf (std)
8. Array smf 
9. Array smf (std)
0. Expected dmf
11. Observed dmf
12. Observed dmf (std)

\* this file does not contain epsilon, but e_scores
instead. You can calculate epsilon from this 
file as follows:

eps = Observed_dmf - Expected_dmf

Final format for 2016 paper
---------------------------
### 11 col
1. Query Strain ID
2. Query allele name
3. Array Strain ID
4. Array allele name
5. Arraytype/Temp
6. Genetic interaction score (ε)
7. P-value
8. Query single mutant fitness (SMF)
9. Array SMF
10. Double mutant fitness
11. Double mutant fitness standard deviation

Final format for the 2010 paper
-------------------------------
### 13 col
1. Query ORF
2. Query gene name
3. Array ORF
4. Array gene name
5. Genetic interaction score (eps)
6. Standard deviation
7. p-value
8. Query single mutant fitness (SMF)
9. Query SMF standard deviation
10. Array SMF
11. Array SMF standard deviation
12. Double mutant fitness
13. Double mutant fitness standard deviation 

