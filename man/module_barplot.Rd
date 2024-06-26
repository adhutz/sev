% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{module_barplot}
\alias{module_barplot}
\title{module_barplot()
Create a bar plot with error bars for the given data}
\usage{
module_barplot(mat, signalome_res)
}
\arguments{
\item{mat}{A matrix containing the data to be plotted}

\item{signalome_res}{A list containing the results of the signalome analysis}
}
\value{
A list of ggplot objects, one for each module in the signalome result
}
\description{
This function takes a matrix and a signalome result, processes the data, and creates a bar plot with error bars.
It returns a list of ggplot objects, one for each module in the signalome result.
}
