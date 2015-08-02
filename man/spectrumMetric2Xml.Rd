\name{spectrumMetric2Xml}

\alias{spectrumMetric2Xml}


\title{Convert a SpecturmMetric into XML for the BSS}

\description{
The \code{spectrumMetric2Xml} function converts a \code{SpectrumMetric} into an
XML structure appropriate for submitting to the MUSTANG Backend Storage System (BSS).
}

\usage{
spectrumMetric2Xml(metricList)
}

\arguments{
  \item{metricList}{a list of \code{SpectrumMetric} objects}
}

%% \details{ }

\value{
A character string with BSS formatted XML is returned.
}

%% \references{ }

\author{
Jonathan Callahan \email{jonathan.s.callahan@gmail.com}
}

%% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{ }

\examples{
  \dontrun{
# Open a connection to IRIS DMC webservices
iris <- new("IrisClient")

# NOTE:  The following trace has 1.728 million points.
# NOTE:  Downloading and calculating PSD may take a while.
starttime <- as.POSIXct("2010-02-27",tz="GMT")
endtime <- as.POSIXct("2010-02-28",tz="GMT")
  
# Get the waveform
st <- getDataselect(iris,"IU","ANMO","00","BHZ",starttime,endtime)

# Make sure we're working with a single snclq
unique_ids <- uniqueIds(st)
if (length(unique_ids) > 1) {
  stop(paste("PSDMetric: Stream has",unique_ids,"unique identifiers"))
}
snclq <- unique_ids[1]

# Calculate and plot the Power Spectral Density
psd <- PSD(st)
  
# Create a Spectrum metric
# TODO:  total revamp of PSD calculations to match Manochehr's needs
phases <- rep(0,length(psd$freq))
m <- new("SpectrumMetric", snclq=snclq, starttime=starttime, endtime=endtime,
          metricName="psd", freqs=psd$freq, amps=psd$spec, phases=psd$phase)

# Show the XML version of the metric
bssXml <- spectrumMetric2Xml(m)
cat(bssXml)
  }
}


