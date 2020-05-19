library(Rsubread)

filepath <- "/Users/danelimjoco/Desktop/GithubRepos/MyHeatmapsWebApp/toy.sam"
counts <- featureCounts(filepath)
counts

deseqdata <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~condition)

