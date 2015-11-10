
PEAKFOREST WEB SERVICE REQUIREMENTS
===================================

 * Get a list of all molecules (IDs only).
 * Get the total number of molecules.
 * Get the names (main name only, i.e. one name per molecule) of molecules. Optional filtering:
    - On a list of molecule IDs.
 * Get a list of all used chromatographic columns. All main information for each column, including column ID or unique name if defined. Optional filtering:
    - On a list of molecule IDs.
 * Find molecules by names. Given a list of names, find for each name the matching molecules (exact case-insensitive match).
 * Get retention times for a molcule ID. Optional filtering:
    - On a list of column.
 * Get the total number of peaks available. Optional filtering:
    - On a list of molecule IDs.
    - On the mode (positive/negative).
 * MS matching on mz/rt. Given a mode (positive/negative), and a mz range, look for all matching peaks. The values to return are: molecule ID, mz, composition and attribution. Optional filtering:
    - On a list of molecule IDs.
    - On a list of attributions.
    - On a range of retention times with a list of columns. In this case the following values are also to be returned: column and matching retention time.

ALGORITHMS
==========

## MZ Annotation

See algorithm \ref{mz}.

\begin{algorithm}
	\caption{MZ Annotation version 1}\label{searchmz}
	\begin{algorithmic}[1]
		\Function{SearchMz}{mzvalues, shift, precision}
			\State $results \gets ()$
			\For{$mz$ in $mzvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$.
				\If{No matched peak}
					\State Append $(mz)$ to $results$.
				\Else
					\ForAll{matching peaks}
						\State Append $(mz, mztheo, molid, composition, attribution)$ to $results$.
					\EndFor
				\EndIf
			\EndFor
			\State \Return results
		\EndFunction
	\end{algorithmic}
\end{algorithm}

## MZ/RT Annotation

See algorithm \ref{searchmzrt}.

\begin{algorithm}
	\caption{MZ/RT Annotation version 1}\label{searchmzrt}
	\begin{algorithmic}[1]
		\Function{SearchMzRt}{mzrtvalues, shift, precision, columns, x, y}
			\State $results \gets ()$
			\For{$(mz, rt)$ in $mzrtvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State $rtinf \gets rt - x - rt^y$
				\State $rtsup \gets rt + x + rt^y$
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$ and $col \in columns$ and $colrt \in [rtinf, rtsup]$.
				\If{No matched peak}
					\State Append $(mz, rt)$ to $results$.
				\Else
					\ForAll{matching peaks}
						\State Append $(mz, rt, mztheo, col, colrt, molid, composition, attribution)$ to $results$.
					\EndFor
				\EndIf
			\EndFor
			\State \Return results
		\EndFunction
	\end{algorithmic}
\end{algorithm}

## MZ Annotation Version 2

In the version 2 of the algorithm, the peak matching is made in two steps:

 1. We search only the peaks whose attribution is among possible precursors, and we keep only the molecule IDs.
 2. We search all peaks, but only in the list of molecules obtains in step 1.

The precursor list is defined as following:

 * In positive mode, $PRECUSOR = \{[(M+H)]^+, [(M+Na)]^+, [(M+K)]^+\}$.
 * In negative mode, $PRECUSOR = \{[(M-H)]^-, [(M+Cl)]^-\}$.

See algorithm \ref{searchmz2}.

\begin{algorithm}
	\caption{MZ Annotation version 2}\label{searchmz2}
	\begin{algorithmic}[1]
		\Function{SearchMz2}{mzvalues, shift, precision}
			\State $molids \gets ()$, a set of unique molecule IDs.
			\For{$mz$ in $mzvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$ and $attribution \in PRECURSORS$.
				\ForAll{matching peaks}
					\State Append $molid$ to the $molids$ set.
				\EndFor
			\EndFor
			\State $results \gets ()$
			\For{$mz$ in $mzvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$ and $molid \in molids$.
				\If{No matched peak}
					\State Append $(mz)$ to $results$.
				\Else
					\ForAll{matching peaks}
						\State Append $(mz, mztheo, molid, composition, attribution)$ to $results$.
					\EndFor
				\EndIf
			\EndFor
			\State \Return results
		\EndFunction
	\end{algorithmic}
\end{algorithm}

## MZ/RT Annotation Version 2

See algorithm \ref{searchmzrt2}.

One issue with this algorithm, is that, because there is no matching on column retention times in the second loop, it returns all columns among the wanted columns, for which we have retention time values. This means that a peak for which we have two of the wanted columns, will give two results, one for each column, even if only one of the columns has been matched in the first loop of the algorithm. We should, inside the $rtmolids$ map, memorize the column matched, and use it in the second loop in order to get only the right match.

\begin{algorithm}
	\caption{MZ/RT Annotation version 2}\label{searchmzrt2}
	\begin{algorithmic}[1]
		\Function{SearchMzRt2}{mzrtvalues, shift, precision, columns, x, y, z}
			\State $rtmolids \gets ()$, a map $rt \rightarrow molids$.
			\For{$(mz, rt)$ in $mzvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State $rtinf \gets rt - x - rt^y$
				\State $rtsup \gets rt + x + rt^y$
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$ and $col \in columns$ and $colrt \in [rtinf, rtsup]$ and $attribution \in PRECURSORS$.
				\State $molids \gets ()$, a set of unique molecule IDs.
				\ForAll{matching peaks}
					\State Append $molid$ to the $molids$ set.
				\EndFor
				\State Append $(rt, molids)$ to the $rtmolids$ map.
			\EndFor
			\State $results \gets ()$
			\For{$(mz, rt)$ in $mzrtvalues$}
				\State $mzinf \gets mz ( 1 + \frac{- shift - precision}{10^6} )$
				\State $mzsup \gets mz ( 1 + \frac{- shift + precision}{10^6} )$
				\State $rtinf \gets rt - x - rt^y$
				\State $rtsup \gets rt + x + rt^y$
				\State $molids \gets ()$, a set of unique molecule IDs.
				\State Append all molecule IDs of $rtmolids$ map whose retention time $\in [rt - z, rt + z]$.
				\State Search for all peaks with $mztheo \in [mzinf, mzsup]$ and $col \in columns$ and $molid \in molids$.
				\If{No matched peak}
					\State Append $(mz, rt)$ to $results$.
				\Else
					\ForAll{matching peaks}
						\State Append $(mz, rt, mztheo, col, colrt, molid, composition, attribution)$ to $results$.
					\EndFor
				\EndIf
			\EndFor
			\State \Return results
		\EndFunction
	\end{algorithmic}
\end{algorithm}
