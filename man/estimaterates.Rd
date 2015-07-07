\name{estimaterates}
\alias{estimaterates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate substitution rate matrix.
}
\description{
Estimate a substitution rate matrix. Users can specify a substitution rate matrix and estimate it. Partition analysis, gamma rate variation, estimating root frequencies of the discrete characters, and specifying clades (or a set of branches) with their own unique rates are also possible. See arguments below.}
\usage{
estimaterates(usertree = NULL, userphyl = NULL, matchtipstodata = FALSE,
              unobserved = NULL, alphabet = NULL, modelmat = NULL, 
              bgtype = "listofnodes", bg = NULL, partition = NULL, 
              ratevar = FALSE, nocat = 4, numhessian = TRUE, 
              rootprob = NULL, rpvec = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{usertree}{
Rooted binary tree of class "phylo". Read in Newick tree using read.tree() from package ape before passing this argument. The branch lengths must be in expected substitutions per site. Trees estimated using MrBayes and BEAST, for instance, yield branch lengths using that scale.
}
  \item{userphyl}{
A matrix of phyletic patterns. Rows represent the discrete character patterns and columns represent the taxa. These data can be numeric or character (but not both). 
}
  \item{matchtipstodata}{
The default is FALSE, which means that the ordering of the taxa in the data matrix must match the internal ordering of tip labels of the tree. Set to TRUE, if the column names of the data matrix, i.e., the taxa names, are all present in the tip labels of the tree provided, in which case, the restriction on ordering is not necessary.
}
  \item{unobserved}{
A matrix of unobserved phyletic patterns, representing possible sampling or acquisition bias. Each row should be a unique phyletic pattern.
}
  \item{alphabet}{
The set of discrete characters used. May be integer only or character only.
}
  \item{modelmat}{
A square matrix denoting the substitution rates to be estimated. The number of rows and columns corresponds to the total number of discete states possible in the data. For example, using matrix(c(1, 2, 2, 1), 2, 2) means that two rates must be estimated corresponding to the entries 1 and 2, respectively. Using matrix matrix(c(1, 2, 0, 1), 2, 2) means that two rates must be estimated as before, but the entry 0 corresponds to a substitution that is not permitted (hence, does not need to be estimated). See examples for more examples.
}
  \item{bgtype}{
Use this to group branches hypothesized to follow the same rates (but differ from other branches). If clade-specific insertion and deletion rates are required to be estimated, use argument option "ancestornodes". If, on the other hand, a group of branches (not in a clade) are hypothesized to follow the same rates, use argument option "listofnodes".}
  \item{bg}{
A vector of nodes should be provided if the "ancestornodes" option was chosen for argument "bgtype". If, on the other hand, "listofnodes" was chosen, a list should be provided with each element of the list being a vector of nodes that limit the branches that follow the same rates. See examples.
}
  \item{partition}{
A list of vectors (of sites) subject to different evolutionary constraints. For example, supplying list(c(1:2500), c(2501:5000)) means that sites 1 through 2500 follow their own substitution rates distinct from sites 2501 through 5000. These sites correspond to the rows of the data supplied in userphyl.
}
  \item{ratevar}{
Default option is FALSE. Specifying "discgamma" implements the discrete gamma approximation model of Yang, 1994. Even if a partition of sites is specified, this option uses the same \eqn{\alpha} parameter for the gamma distribution over all partitions. Specifying "partitionspecificgamma" implements the discrete gamma approximation model in each partition separately.
The number of categories can be specified using the "nocat" argument. See examples.
}
  \item{nocat}{
The number of categories for the discrete gamma approximation. 
}
  \item{numhessian}{
Set to FALSE if standard errors are not required. This speeds up the algorithm. Default is TRUE.
}
  \item{rootprob}{
Four options are available: "equal", "stationary", "maxlik", and "user". 
Option "equal" means that all the discrete characters are given equal weight at the root. 
Using "stationary" means that the discrete characters are weighted at the root by the stationary frequencies implied by the substitution rate matrix. Note that these can differ based on the partition. In the case that the "bgtype" argument is also provided, an average of the stationary frequencies of all branch groupings is used.
If "maxlik" is supplied, then the root frequencies are also estimated. These do not differ based on the partition.
If "user" is supplied here, a vector of root frequency parameters can be provided to argument "rpvec".
}
  \item{rpvec}{
If option "user" is specified for argument "rootprob", supply a vector of the same length as that provided to argument "alphabet", representing the root frequency parameters.
}
  \item{\dots}{
Passing other arguments to the optimization algorithm "nlminb". For example, control = list(trace = 5) will print progress at every 5th iteration.}
}
\details{
The discrete character patterns should be ordered by the tips of the tree.
}
\value{
All arguments used while calling the estimaterates function are attached in a list. The following components are also returned in the same list:
\item{call        }{Function call used.}
\item{conv        }{A vector of convergence indicators for the model run. 0 denotes successful convergence.}
\item{time        }{Time taken in seconds.}
\item{tree        }{The phylogenetic tree used.}
\item{bg          }{List of group of nodes that capture branches that follow unique substitution rates.}
\item{results     }{List of results including the parameter estimates for each unique entry of modelmat (excluding zeros), standard errors, number of parameters fit, and AIC and BIC values. Furthermore, details from the optimization routine applied are also available.}
\item{data_red    }{Unique phyletic patterns observed.}
\item{w           }{Number of times each gene phyletic pattern was observed.}
}

\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@mcmaster.ca}
}
\note{
See the vignette for an in-depth look at some examples.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#library(markophylo)
##############
#data(simdata1) #example data for a 2-state continuous time markov chain model.
#Now, plot example tree.
#ape::plot.phylo(simdata1$tree, edge.width = 2, show.tip.label = FALSE, no.margin = TRUE) 
#ape::nodelabels(frame = "circle", cex = 0.7)
#ape::tiplabels(frame = "circle", cex = 0.7)
#print(simdata1$Q) #substitution matrix used to simulate example data
#print(table(simdata1$data)) #states and frequencies
#model1 <- estimaterates(usertree = simdata1$tree, userphyl = simdata1$data, 
#                        alphabet = c(1, 2), rootprob = "equal", 
#                        modelmat = matrix(c(NA, 1, 2, NA), 2, 2))
#print(model1)
####
#If the data is known to contain sampling bias such that certain phyletic
#patterns are not observed, then these unobserved data can be corrected for
#easily. First, let's create a filtered version of the data following which
#a correction will be applied within the function.
#filterall1 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
#                            function(x) isTRUE(all.equal(as.vector(x), c(1, 1, 1, 1)))))
#filterall2 <- which(apply(simdata1$data, MARGIN = 1, FUN = 
#                            function(x) isTRUE(all.equal(as.vector(x), c(2, 2, 2, 2)))))
#filteredsimdata1 <- simdata1$data[-c(filterall1, filterall2), ]
#model1_f_corrected <- estimaterates(usertree = simdata1$tree, userphyl = filteredsimdata1, 
#                          unobserved = matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 2,  byrow = TRUE), 
#                          alphabet = c(1, 2), rootprob = "equal", 
#                          modelmat = matrix(c(NA, 1, 2, NA), 2, 2))
#print(model1_f_corrected)
##############
#data(simdata2)
#print(simdata2$Q)
#While simulating the data found in simdata2, the clade with node 7 as its
#most recent common ancestor (MRCA) was constrained to have twice the 
#substitution rates as the rest of the branches in the tree.
#print(table(simdata2$data))
#model2 <- estimaterates(usertree = simdata2$tree, userphyl = simdata2$data, 
#                        alphabet = c(1, 2), bgtype = "ancestornodes", bg = c(7),
#                        rootprob = "equal", modelmat = matrix(c(NA, 1, 2, NA), 2, 2))
#print(model2)
#plottree(model2, colors=c("blue", "darkgreen"), edge.width = 2, show.tip.label = FALSE, 
#         no.margin = TRUE)
#ape::nodelabels(frame = "circle", cex = 0.7)
#ape::tiplabels(frame = "circle", cex = 0.7)
##############
#Nucleotide data was simulated such that the first half of sites followed
#substitution rates different from the other half of sites. Data was simulated
#in the two partitions with rates 0.33 and 0.99.
#data(simdata3)
#print(dim(simdata3$data))
#print(table(simdata3$data))
#model3 <- estimaterates(usertree = simdata3$tree, userphyl = simdata3$data, 
#                        alphabet = c("a", "c", "g", "t"), rootprob = "equal", 
#                        partition = list(c(1:2500), c(2501:5000)), 
#                        modelmat = matrix(c(NA, 1, 1, 1, 1, NA, 1, 1, 
#                                            1, 1, NA, 1, 1, 1, 1, NA), 4, 4))
#print(model3)
#
#More examples in the vignette.
}