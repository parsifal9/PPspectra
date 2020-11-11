We  consider the question of analysing a point process (a collection of "times" in the frequency domain.
This raises the question of how, if we model a point process as series of spikes at particular times, we transform the
time domain to the frequence domain. We mention here that:
\begin{itemize}
\item there are many excellent packages in R \cite[]{r.2013}, for the time domain analysis of point processes. There
  appears to be a dearth of software for the frequence analysis of point processes;
\item the choice to model the data as a "spike train" may be more natural in some settings than others.  The
  term "spike train" comes from the neuroscience literature.  \cite[]{Eden.2019,Pouzat.and.Chaffiol.2009} model the firing
  times of neurons as a point process that is then analysed as a ``spike train'' where each event is an electrical
  pulse, a signal of zero time duration with some particular mathematical properties. 
\end{itemize}

The frequency domain has several advantages in that subtle structure can be detected which may be difficult to observe
in the time domain. Secondly, the time domain estimators can have problems with sensitivity to weak non-stationarity.
\cite{Jarvis.and.Mitra.2001} is a key reference in this area. As they comment, most of the literature is targeted at
either spectral analysis of continuous processes or at the analysis of point processes but in the time domain.

We could find no package in R that explicitly offered the calculation of an energy spectra for a point process.  The
Chronux Matlab package \cite[]{Mitra.and.Bokil.2007,chronux}, offered this so we have taken a small number of functions
and coded them in R. We have kept the same names ( \texttt{mtspectrumpt} and \texttt{coherencypt}) for uses familiar
with Chronux, but the calling parameters may have changed.

```
library(devtools)
devtools::install_github("parsifal9/PPspectra", build_vignettes = TRUE)
```
