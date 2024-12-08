% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{corr_plot}
\alias{corr_plot}
\title{Create a correlation plot with different conditions and labels}
\usage{
corr_plot(
  df,
  x_column,
  y_column,
  label_names,
  gene_list = c(""),
  top_genes = 20,
  max.overlaps = Inf
)
}
\arguments{
\item{df}{A data frame with columns for x, y, and label names.}

\item{x_column}{Character string specifying the column name for the x-axis values.}

\item{y_column}{Character string specifying the column name for the y-axis values.}

\item{label_names}{Character string specifying the column name for the labels.}

\item{gene_list}{A character vector of gene names to be highlighted. Default is an empty vector.}

\item{top_genes}{Integer specifying the number of top differing proteins to consider. Default is  20.}

\item{max.overlaps}{Integer specifying the maximum number of overlapping labels. Default is Inf.}
}
\value{
A list containing the plot and subsets of the data used for each part of the plot.
}
\description{
This function takes a dataframe and creates a scatter plot with correlation
analysis, highlighting points based on certain criteria such as being in a
gene list or having a rank among the top N differences.
}
