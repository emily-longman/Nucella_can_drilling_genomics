# This is known to work at least under the following conditions:
# - Oxford Nanopore Ultra-Long reads.
# - Human genomes.
# - Guppy 5 base caller or newer with "super" accuracy.
# - Coverage 40x to 80x. If you have more coverage, 
#   adjust --Reads.minReadLength or --Reads.desiredCoverage
#   to bring coverage down to this range.

# Under the above conditions, this should give an assembly with
# N50 in the tens of Mb.


[Reads]
minReadLength = 5000
noCache = True

[Kmers]
k = 14

[MinHash]
minBucketSize = 5
maxBucketSize = 30
minFrequency = 5

[Align]
alignMethod = 3
downsamplingFactor = 0.05
matchScore = 6
sameChannelReadAlignment.suppressDeltaThreshold = 30

# The following Align parameters are set to very permissive values to allow the majority of alignments
# to be assessed during the initial stage of automatic alignment parameter selection
# (ReadGraph.creationMethod 2).
maxSkip = 100
maxDrift = 100
maxTrim = 100
minAlignedMarkerCount = 10
minAlignedFraction = 0.1

[ReadGraph]
# This uses the observed distribution of alignment statistics to choose thresholds for
# maxSkip, maxDrift, maxTrim, minAlignedMarkerCount, and minAlignedFraction
creationMethod = 2

[MarkerGraph]
simplifyMaxLength = 10,100,1000,10000,100000
crossEdgeCoverageThreshold = 3

# Adaptive estimation of coverage threshold to generate marker graph vertices.
minCoverage = 0

[Assembly]
consensusCaller = Bayesian:guppy-5.0.7-b
detangleMethod = 2