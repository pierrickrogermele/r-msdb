% RMsDb
% Pierrick Roger Mele
% March 3, 2015

# Environment

## Used technology

 * R language.

 * Direct connection to PostgreSQL database.

## Class diagram

<center>![Class diagram](diagrams/msdb-class-diagram.png)</center>

# Input format

 * TSV (Tabulation separated values) file for the R script.
 * Excel file for the website.

<center>
<!-- The empty column is here to insert some space between the two columns. -->

MZ                      RT
--------------  -----   ---------------
75.02080998             49.38210915
75.05547146             0.658528069
75.08059797             1743.94267
76.03942694             51.23158899
...                     ...

</center>

# Output format

 * TSV (Tabulation separated values) file for the R script.
 * Excel file for the website.

<center>
<!-- The empty column is here to insert some space between the two columns. -->

MZ	                    RT	                MZTHEO	    COL	        COLRT	ID	COMPOSITION	    ATTRIBUTION	      
--------------  -----   ---------------     ----------  ----------- -----   --- ------------    ------------------
76.03942694	            51.23158899	        76.039305	"UPLC (C8)"	47.4	156	"C2 H6 N O2"	"[(M+H)-(C2H5N)]+" 
76.07584477	            50.51249853	        76.07569	"UPLC (C8)"	50.4	471	"C3 H10 N O"	"[(M+H)-(CO2)]+"   
76.07584477	            50.51249853	        76.07569	"UPLC (C8)"	67.8	233	"C3 H10 N O"	"[(M+H)]+"	       
76.07593168	            0.149308136	        NA	        NA	        NA	    NA	NA	            NA	               
...                     ...

</center>

# Annot 1 Without RT

## Parameters:

 * `shift`, default value is 0.
 * `precision`, default value is 5.

## Range of search:

 * $minmz = mz (1 + \frac{- shift - precision}{10^6})$.
 * $maxmz = mz (1 + \frac{- shift + precision}{10^6})$.

# Annot 1 With RT

## Parameters:

 * `col`, a list of columns.
 * `x`, default value is 5.
 * `y`, default value is 0.8.

## Range of search:

 * $minrt = rt - x - rt^y$
 * $maxrt = rt + x + rt^y$

# Annot 2 Without RT

## Parameters:

 * `shift`, default value is 1.
 * `precision`, default value is 10.

## Search

This second version of the algorithm works in several steps:

 1. It looks only at peaks whose attribution is a precusor: [(M+H)]+, [(M+Na)]+, [(M+Cl)]-, ...
 2. It constructs a list of matched molecules, using the results of step 1.
 3. It looks at all peaks inside the molecule list, using the same range equations as in version 1.

# Annot 2 With RT

## Parameters:

 * `col`, a list of columns.
 * `x`, default value is 5.
 * `y`, default value is 0.8.
 * `z`, default value is 5.

## Search

This second version of the algorithm works in two steps:

 1. It looks only at peaks whose attribution is a precusor ([(M+H)]+, [(M+Na)]+, [(M+Cl)]-, etc), using retention time equations with `x` and `y`.
 2. It constructs a map of retention time / molecules, from the results of step 1.
 3. It looks at all peaks inside a molecule list built from elements of the map whose retention time is between $rt - z$ and $rt + z$.
