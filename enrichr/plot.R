options(stringsAsFactors = FALSE)
library(methods)
library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
	cat("Usage: plot.R <enrichr result file> <output file> <query size> <resolution> <width> <height>\n", file = stderr())
	q(save = 'no')
}

erfile  = args[1]
outfile = args[2]
qsize   = as.numeric(args[3])
res     = as.numeric(args[4])
width   = as.numeric(args[5])
height  = as.numeric(args[6])

erdata = read.table(erfile, header = TRUE, row.names = NULL, check.names = FALSE, sep = "\t", quote = '"')
n = nrow(erdata)
if (n == 0) {
	cat("No terms are enriched, nothing to plot.", file = stderr())
	q(save = "no")
}
terms = sapply(c(erdata$Term), function(x) {
    if (nchar(x) > 45) {
        occurs = unlist(gregexpr("[_(]", x))
        lastoc = occurs[length(occurs)]
        suffix = substring(x, lastoc)
        prefix = substr(x, 1, lastoc - 1)
        paste(substr(x, 1, 45 - (nchar(x) - lastoc + 1)), '...', suffix)
    } else x
})
ovdata = t(as.data.frame(strsplit(as.character(erdata$Overlap), '/', fixed=T)))
pdata = data.frame(
    Term = factor(terms, levels=rev(unique(terms))),
    LogPval = -log10(erdata$AdjPval),
    OverlapN = as.numeric(ovdata[, 1]),
    TermN = sprintf("n = %-3s", ovdata[, 2])
)

maxOvN = max(pdata$OverlapN)
maxLogPval = max(pdata$LogPval, 1)
pdata$LogPval = pdata$LogPval/maxLogPval*maxOvN
pcut = -log10(.05) * maxOvN / maxLogPval
pvalticks = ceiling(max(maxLogPval, -log10(.05)))
if (pvalticks >=6) {
	ptbreaks = 0:pvalticks
	ptlabels = 0:pvalticks
} else if ( pvalticks >= 3 ) {
    pvalticks = 0:pvalticks
	ptbreaks  = c(pvalticks, -log10(.05))
	ptlabels  = c(pvalticks, 'p=0.05')
} else {
    # remove 1.5, cuz it's too close to 1.3 (-log10(.05))
    pvalticks =  c(0, .5, 1, seq(2, pvalticks, by = .5))
	ptbreaks  = c(pvalticks, -log10(.05))
	ptlabels  = c(pvalticks, 'p=0.05')
}

library(ggplot2)
p = ggplot(pdata) +
	geom_col(aes(x = Term, y = OverlapN), fill = "blue", color = '#4286f4', alpha = .1) +
    geom_hline(yintercept = pcut, color = "red", alpha = .3, linetype = 'dashed') +
	geom_line(aes(x = Term, y = LogPval, group = 1), color = "red", size = .5)  +
	geom_point(aes(x = Term, y = LogPval, group = 1), fill = "white", color = "red3", size = 2, shape = 21) +
    geom_text(aes(x = Term, y = 0, label = TermN, group = 1), hjust = -.2, size = 3) +
	scale_y_continuous(
        name = sprintf("# overlapping genes (query size: %s)", qsize),
		sec.axis = sec_axis(
            ~ . * maxLogPval / maxOvN,
            name = "-log10(adj.pval)",
            breaks = ptbreaks,
            labels = ptlabels
        )
	) +
	coord_flip() + theme_bw() +
	theme(
		panel.background = element_rect(fill = alpha("#EFEFEF", .4)),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line.x.top = element_line(color = "red3", size = .6),
		axis.ticks.x.top = element_line(color = "red3"),
		axis.line.x.bottom = element_line(color = "blue3", size = .6),
		axis.ticks.x.bottom = element_line(color = "blue3"),
		axis.title.y.left = element_blank()
	)
# scale height
# n = 20, width = 2000, height = 1900
height = as.integer(1900.0/2000*width*n/20)
png(outfile, res = res, width = width, height = height)
print(p)
dev.off()

