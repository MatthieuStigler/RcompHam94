\subsection{Fractional Integration}
This example uses package \lib{fracdiff} to compute the exponent of fractional integration as described on pp 448-449. We use the function \fun{fdGPH}:
\begin{Scode}
library(fracdiff)
args(fdGPH)
\end{Scode}

Applied on US GDP and Treasury Yields data:

\begin{Scode}{fracdiff}
data(gnptbill,package = "RcompHam94")
print( fdGPH(gnptbill$GNP) )
print( fdGPH(gnptbill$TBILL) )
\end{Scode}
