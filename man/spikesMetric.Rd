\name{spikesMetric}

\alias{spikesMetric}


\title{Find spikes using a rolling Hampel filter}

\description{
The spikesMetric() function determines the number of spikes in a seismic \code{Stream}.
}

\usage{
spikesMetric(st, windowSize, thresholdMin, selectivity)
}

\arguments{
  \item{st}{ a \code{Stream} object containing a seismic signal}
  \item{windowSize}{ The window size to roll over (default=\code{41}) }
  \item{thresholdMin}{ Minimum value for outlier detection (default=\code{6.0}) }
  \item{selectivity}{ Numeric factor [0-1] used in determining outliers (default=\code{0.4}) }
}

\details{
This function uses the output of the \code{findOutliers} function in the \pkg{seismicRoll}
package to calculate the number of 'spikes' containing outliers.

The \code{thresholdMin} level is similar to a sigma value for normally distributed data.
Hampel filter values above 6.0 indicate a data value that is extremely unlikely
to be part of a normal distribution  (~ 1/500 million) and therefore very likely to be an outlier. By
choosing a relatively large value for \code{thresholdMin} we make it less likely that we
will generate false positives.

The \code{selectivity} is a value between 0 and 1 and is used to generate an appropriate 
threshold for outlier detection based on the statistics of the incoming data. A lower value
for \code{selectivity} will result in more outliers while a value closer to 1.0 will result in 
fewer.

The total count of spikes reflects the number of outlier data points that are separated by at least 
one non-outlier data point. Each individual spike may contain more than one data point.

}

\value{
A list of \code{SingleValueMetric} objects is returned.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

\note{
The \code{thresholdMin} parameter is sensitive to the data sampling rate. The default
value of 6.0 seems to work well with sampling rates of 10 Hz or higher ('B..' or 'V..' channels).  For
'L..' channels with a sampling rate of 1 Hz it seems that \code{thresholdMin=12.0} is more appropriate.

\strong{More testing of spiky signals at different resolutions is needed.}

%%  ~~further notes~~
See the \pkg{seismicRoll} package for documentation on the findOutliers() function.
}

%% \seealso{ %% ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# Get the waveform
starttime <- as.POSIXct("2013-01-03 15:00:00", tz="GMT")
endtime <- starttime + 3600 * 3  
st <- getDataselect(iris,"IU","RAO","10","BHZ",starttime,endtime)

# Calculate the gaps metrics and show the results
metricList <- spikesMetric(st)
dummy <- lapply(metricList, show)
}

\keyword{metrics}
