\name{crossCorrelationMetrics}

\alias{crossCorrelationMetrics}


\title{Correlation between channels}

\description{
The crossCorrelationMetrics() function calculates the \code{peak_correlation}
and \code{peak_lag} associated with two streams of seismic data.
}

\usage{
crossCorrelationMetrics(st1, st2, maxLagSecs=10, filter=signal::butter(2,0.2))
}

\arguments{
  \item{st1}{a \code{Stream} object containing a seismic signal}
  \item{st2}{a \code{Stream} object containing a seismic signal}
  \item{maxLagSecs}{maximum number of seconds of lag to use}
  \item{filter}{a \pkg{signal} package filter to be applied before cross-correlating}
}

\details{
%% max preformatted length ---------------------------------------------------
Details of the algorithm are as follows:
\itemize{
  \item{Both signals are demeaned and detrended}
  \item{If one signal has a higher sampling rate, it is decimated to the lower sampling rate using an IIR filter see(\code{signal::decimate}).}
  \item{Both signals are filtered, by default with a Butterworth low pass filter see(\code{signal::filter}).}
  \item{Signals are cross-correlated using the ccf() function.}
}

The maximum absolute correlation is saved as \code{peak_correlation} while the lag at peak correlation is saved as \code{peak_lag}.

\bold{Note:} For cross-correlation, seismic signals must not have any gaps -- they must be contained in a single \code{Trace} object.
}

\value{
A list with two \code{SingleValueMetric} objects is returned.
The metric names are \code{peak_correlation} and \code{peak_lag}.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com} (R code),
Mary Templeton \email{met@iris.washington.edu} (algorithm)
}

\note{
The metricList generated for this two-channel metric will have an additional \code{sncl2} attribute identifying the SNCL in \code{st2}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{ }

\examples{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# Get the same signal, shifted by 3 seconds
starttime <- as.POSIXct("2013-11-12 07:09:45",tz="GMT")
endtime <- starttime + 600
st1 <- getSNCL(iris,"NM.SLM.00.BHZ",starttime,endtime)
st2 <- getSNCL(iris,"NM.SLM.00.BHZ",starttime+3,endtime+3)

# Cross-correlate
crossCorrelationMetrics(st1,st2)

#[[1]]
# SingleValueMetric 
#  metric:        peak_correlation
#  snclq:         NM.SLM.00.BHZ.M
#  starttime:     2013-11-12 07:09:45
#  endtime:       2013-11-12 07:19:45
#  value:         0.999
#  snclq2:    NM.SLM.00.BHZ.M
#
#[[2]]
#SingleValueMetric 
#  metric:        peak_lag
#  snclq:         NM.SLM.00.BHZ.M
#  starttime:     2013-11-12 07:09:45
#  endtime:       2013-11-12 07:19:45
#  value:         3.000
#  snclq2:    NM.SLM.00.BHZ.M
}

\keyword{metrics}
