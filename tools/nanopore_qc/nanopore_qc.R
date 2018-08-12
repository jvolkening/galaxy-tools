#!/usr/bin/Rscript

# MinionQC version 1.0
# Copyright (C) 2017 Robert Lanfear
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# supress warnings
options(warn=-1)

library(ggplot2)
library(plyr)
library(reshape2)
library(readr)
library(yaml)
suppressMessages(library(scales))
#library(parallel)
library(futile.logger)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))


# option parsing #

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-i", "--input"), 
                     type = "character",
                     dest = 'input.file',
                     help="Input file or directory (required). Either a full path to a sequence_summary.txt file, or a full path to a directory containing one or more such files. In the latter case the directory is searched recursively."
                     )

parser <- add_option(parser, 
                     opt_str = c("-o", "--outputdirectory"), 
                     type = "character",
                     dest = 'output.dir',
                     help="Output directory (required). If a single sequencing_summary.txt file is passed as input, then the output directory will contain just the plots associated with that file. If a directory containing more than one sequencing_summary.txt files is passed as input, then the plots will be put into sub-directories that have the same names as the parent directories of each sequencing_summary.txt file"
                     )

parser <- add_option(parser, 
                     opt_str = c("-q", "--qscore_cutoff"), 
                     type="double", 
                     default=7.0,
                     dest = 'q',
                     help="The cutoff value for the mean Q score of a read (default 7). Used to create separate plots for reads above and below this threshold"
                     )

parser <- add_option(parser, 
                     opt_str = c("-d", "--discard_failed"), 
                     type="logical", 
                     default=FALSE,
                     dest = 'filt.failed',
                     help="Discard reads that failed Albacore filtering"
                     )

opt = parse_args(parser)

input.file  = opt$input.file
output.dir  = opt$output.dir
filt.failed = opt$filt.failed
q           = opt$q

# this is how we label the reads at least as good as q
q_title = paste("Q>=", q, sep="")


# build the map for R9.5
p1 = data.frame(channel=33:64, row=rep(1:4, each=8), col=rep(1:8, 4))
p2 = data.frame(channel=481:512, row=rep(5:8, each=8), col=rep(1:8, 4))
p3 = data.frame(channel=417:448, row=rep(9:12, each=8), col=rep(1:8, 4))
p4 = data.frame(channel=353:384, row=rep(13:16, each=8), col=rep(1:8, 4))
p5 = data.frame(channel=289:320, row=rep(17:20, each=8), col=rep(1:8, 4))
p6 = data.frame(channel=225:256, row=rep(21:24, each=8), col=rep(1:8, 4))
p7 = data.frame(channel=161:192, row=rep(25:28, each=8), col=rep(1:8, 4))
p8 = data.frame(channel=97:128, row=rep(29:32, each=8), col=rep(1:8, 4))

q1 = data.frame(channel=1:32, row=rep(1:4, each=8), col=rep(16:9, 4))
q2 = data.frame(channel=449:480, row=rep(5:8, each=8), col=rep(16:9, 4))
q3 = data.frame(channel=385:416, row=rep(9:12, each=8), col=rep(16:9, 4))
q4 = data.frame(channel=321:352, row=rep(13:16, each=8), col=rep(16:9, 4))
q5 = data.frame(channel=257:288, row=rep(17:20, each=8), col=rep(16:9, 4))
q6 = data.frame(channel=193:224, row=rep(21:24, each=8), col=rep(16:9, 4))
q7 = data.frame(channel=129:160, row=rep(25:28, each=8), col=rep(16:9, 4))
q8 = data.frame(channel=65:96, row=rep(29:32, each=8), col=rep(16:9, 4))

map = rbind(p1, p2, p3, p4, p5, p6, p7, p8, q1, q2, q3, q4, q5, q6, q7, q8)

add_cols <- function(d, min.q){
    # take a sequencing sumamry file (d), and a minimum Q value you are interested in (min.q)
    # return the same data frame with the following columns added 
        # cumulative.bases
        # hour of run
        # reads.per.hour
    
    d = subset(d, mean_qscore_template >= min.q)
    
    if(nrow(d)==0){ 
        flog.error(paste("There are no reads with a mean Q score higher than your cutoff of ", min.q, ". Please choose a lower cutoff and try again.", sep = ""))
        quit()
    }
    
    d = merge(d, map, by="channel")
    d = d[with(d, order(start_time)), ] # sort by read length
    d$cumulative.bases.time = cumsum(as.numeric(d$start_time))
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length
    d$cumulative.bases = cumsum(as.numeric(d$sequence_length_template))
    d$hour = d$start_time %/% 3600
    
    # add the reads generated for each hour
    reads.per.hour = as.data.frame(table(d$hour))
    names(reads.per.hour) = c("hour", "reads_per_hour")
    reads.per.hour$hour = as.numeric(as.character(reads.per.hour$hour))
    d = merge(d, reads.per.hour, by = c("hour"))    
    return(d)
}


load_summary <- function(filepath, min.q){
    # load a sequencing summary and add some info
    # min.q is a vector of length 2 defining 2 levels of min.q to have
    # by default the lowest value is -Inf, i.e. includes all reads. The 
    # other value in min.q is set by the user at the command line
    d = read_tsv(filepath, col_types = cols_only(channel = 'i', 
                                                passes_filtering = 'c',
                                                num_events_template = 'i', 
                                                sequence_length_template = 'i', 
                                                mean_qscore_template = 'n',
                                                sequence_length_2d = 'i',
                                                mean_qscore_2d = 'n',
                                                start_time = 'n'))
    
    if("sequence_length_2d" %in% names(d)){
        # it's a 1D2 or 2D run
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_2d))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_2d))
        d$num_events_template = NA
        d$start_time = as.numeric(as.character(d$start_time))
        
    }else{
        d$sequence_length_template = as.numeric(as.character(d$sequence_length_template))
        d$mean_qscore_template = as.numeric(as.character(d$mean_qscore_template))
        d$num_events_template = as.numeric(as.character(d$num_events_template))
        d$start_time = as.numeric(as.character(d$start_time))
    }

    # ignore 0-length reads
    d <- d[d$sequence_length_template > 0,]
    # ignore reads failing filtering
    if (filt.failed) {
        d <- d[d$passes_filtering == 'True',]
    }
        
    d$events_per_base = d$num_events_template/d$sequence_length_template

    flowcell = basename(dirname(filepath))
    
    # add columns for all the reads
    d1 = add_cols(d, min.q[1])
    d1$Q_cutoff = "All reads"
    
    # add columns for just the reads that pass the user Q threshold
    d2 = add_cols(d, min.q[2])
    d2$Q_cutoff = q_title
    
    # bind those two together into one data frame
    d = as.data.frame(rbindlist(list(d1, d2)))

    # name the flowcell (useful for analyses with >1 flowcell)
    d$flowcell = flowcell
    
    # make sure this is a factor
    d$Q_cutoff = as.factor(d$Q_cutoff)

    
    keep = c("hour","start_time", "channel", "sequence_length_template", "mean_qscore_template", "row", "col", "cumulative.bases", "cumulative.bases.time", "reads_per_hour", "Q_cutoff", "flowcell", "events_per_base")
    d = d[keep]

    d$start_bin = cut(d$start_time, 9,labels=c(1:9))
        
    return(d)
}

reads.gt <- function(d, len){
    # return the number of reads in data frame d
    # that are at least as long as length len
    return(length(which(d$sequence_length_template>=len)))
}

bases.gt <- function(d, len){
    # return the number of bases contained in reads from 
    # data frame d
    # that are at least as long as length len
    reads = subset(d, sequence_length_template >= len)
    return(sum(as.numeric(reads$sequence_length_template)))
}

log10_minor_break = function (...){
    # function to add minor breaks to a log10 graph
    # hat-tip: https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
    function(x) {
        minx         = floor(min(log10(x), na.rm=T))-1;
        maxx         = ceiling(max(log10(x), na.rm=T))+1;
        n_major      = maxx-minx+1;
        major_breaks = seq(minx, maxx, by=1)
        minor_breaks = 
            rep(log10(seq(1, 9, by=1)), times = n_major)+
            rep(major_breaks, each = 9)
        return(10^(minor_breaks))
    }
}

binSearch <- function(min, max, df, t = 100000) {
    # binary search algorithm, thanks to https://stackoverflow.com/questions/46292438/optimising-a-calculation-on-every-cumulative-subset-of-a-vector-in-r/46303384#46303384
    # the aim is to return the number of reads in a dataset (df)
    # that comprise the largest subset of reads with an N50 of t
    # we use this to calculte the number of 'ultra long' reads
    # which are defined as those with N50 > 100KB
    mid = floor(mean(c(min, max)))
    if (mid == min) {
        if (df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[min]/2))] < t) {
            return(min - 1)
        } else {
            return(max - 1)
        }
    }
    
    n = df$sequence_length_template[min(which(df$cumulative.bases>df$cumulative.bases[mid]/2))]
    if (n >= t) {
        return(binSearch(mid, max, df))
    } else {
        return(binSearch(min, mid, df))
    }
}


summary.stats <- function(d, Q_cutoff="All reads"){
    # Calculate summary stats for a single value of min.q
    
    rows = which(as.character(d$Q_cutoff)==Q_cutoff)
    d = d[rows,]
    d = d[with(d, order(-sequence_length_template)), ] # sort by read length, just in case
    
    total.bases = sum(as.numeric(d$sequence_length_template))
    total.reads = nrow(d)
    N50.length = d$sequence_length_template[min(which(d$cumulative.bases > (total.bases/2)))]
    mean.length = round(mean(as.numeric(d$sequence_length_template)), digits = 1)
    median.length = round(median(as.numeric(d$sequence_length_template)), digits = 1)
    max.length = max(as.numeric(d$sequence_length_template))
    mean.q = round(mean(d$mean_qscore_template), digits = 1)
    median.q = round(median(d$mean_qscore_template), digits = 1)
    
    #calculate ultra-long reads and bases (max amount of data with N50>100KB)
    ultra.reads = binSearch(1, nrow(d), d, t = 100000)    
    if(ultra.reads>=1){
        ultra.gigabases = sum(as.numeric(d$sequence_length_template[1:ultra.reads]))/1000000000
    }else{
        ultra.gigabases = 0
    }
        
    reads = list(
                reads.gt(d, 10000), 
                reads.gt(d, 20000), 
                reads.gt(d, 50000),
                reads.gt(d, 100000),
                reads.gt(d, 200000),
                reads.gt(d, 500000),
                reads.gt(d, 1000000),
                ultra.reads)
    names(reads) = c(">10kb", ">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")

    bases = list(
                bases.gt(d, 10000)/1000000000, 
                bases.gt(d, 20000)/1000000000, 
                bases.gt(d, 50000)/1000000000,
                bases.gt(d, 100000)/1000000000,
                bases.gt(d, 200000)/1000000000,
                bases.gt(d, 500000)/1000000000,
                bases.gt(d, 1000000)/1000000000,
                ultra.gigabases)
    names(bases) = c(">10kb", ">20kb", ">50kb", ">100kb", ">200kb", ">500kb", ">1m", "ultralong")
    
    return(list('total.gigabases' = total.bases/1000000000,
                'total.reads' = total.reads,
                'N50.length' = N50.length, 
                'mean.length' = mean.length, 
                'median.length' = median.length,
                'max.length' = max.length,
                'mean.q' = mean.q,
                'median.q' = median.q,
                'reads' = reads,
                'gigabases' = bases
                ))
}

channel.summary <- function(d){
    # calculate summaries of what happened in each of the channels 
    # of a flowcell
    
    a = ddply(d, .(channel), 
              summarize, 
              total.bases = sum(sequence_length_template), 
              total.reads = sum(which(sequence_length_template>=0)), 
              mean.read.length = mean(sequence_length_template), 
              median.read.length = median(sequence_length_template))
    b = melt(a, id.vars = c("channel"))
    return(b)    
}


single.flowcell <- function(input.file, output.dir, q=8){
    # wrapper function to analyse data from a single flowcell
    # input.file is a sequencing_summary.txt file from a 1D run
    # output.dir is the output directory into which to write results
    # q is the cutoff used for Q values, set by the user
    flog.info(paste("Loading input file:", input.file))
    d = load_summary(input.file, min.q=c(-Inf, q))

    flowcell = unique(d$flowcell)

    flog.info(paste(sep = "", flowcell, ": creating output directory:", output.dir))
    dir.create(output.dir)
    out.txt = file.path(output.dir, "summary.yaml")
    
    flog.info(paste(sep = "", flowcell, ": summarising input file for flowcell"))
    all.reads.summary = summary.stats(d, Q_cutoff = "All reads")
    q10.reads.summary = summary.stats(d, Q_cutoff = q_title)
    
    summary = list("input file" = input.file,
                   "All reads" = all.reads.summary,
                   cutoff = q10.reads.summary,
                   "notes" = 'ultralong reads refers to the largest set of reads with N50>100KB')
    
    names(summary)[3] = q_title
    
    write(as.yaml(summary), out.txt)
    
    muxes = seq(from = 0, to = max(d$hour), by = 8)
    
    # make plots
    flog.info(paste(sep = "", flowcell, ": plotting length histogram"))
    len.lims <- quantile(d$sequence_length_template,c(0.01,0.99))
    len.lims[1] <- max(1,len.lims[1])
    
    p1 = ggplot(d, aes(x = sequence_length_template)) + 
        geom_histogram(bins = 200, fill="steelblue") + 
        scale_x_log10(breaks=round(10^pretty(log10(len.lims),n=5),0),limits=len.lims) + 
        facet_wrap(~Q_cutoff, ncol = 1, scales = "free_y") + 
        #theme(text = element_text(size = 15)) +
        theme_linedraw() +
        xlab("Read length") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "length_histogram.png"), width = 7, height = 4, dpi=300, plot = p1)
    ggsave(filename = file.path(output.dir, "length_histogram.screen.png"), width = 7, height = 4, dpi=90, plot = p1)

    flog.info(paste(sep = "", flowcell, ": plotting mean Q score histogram"))
    q.lims <- quantile(d$mean_qscore_template,c(0.005,0.995))
    p2 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = mean_qscore_template)) + 
        geom_histogram(bins = 200, fill="steelblue") + 
        xlim(q.lims) +
        geom_vline(xintercept=q) +
        #theme(text = element_text(size = 15)) +
        theme_linedraw() +
        xlab("Mean Q score of read") +
        ylab("Number of reads")
    ggsave(filename = file.path(output.dir, "q_histogram.png"), width = 7, height = 4, dpi=300, plot = p2)
    ggsave(filename = file.path(output.dir, "q_histogram.screen.png"), width = 7, height = 4, dpi=90, plot = p2)

   
    flog.info(paste(sep = "", flowcell, ": plotting flowcell overview"))

    df <- data.frame(row=integer(),col=integer(),bin=integer(), quality=numeric())

    # swap rows and columns
    for (b in c(1,5,9)) {
        m <- matrix(data=rep(0,512),nrow=16,ncol=32)
        for (r in 1:32) {
        for (c in 1:16) {
            sub <- subset(d, Q_cutoff=="All reads" & col == c & row == r & start_bin == b)
            Q <-  median(sub$mean_qscore_template)
            df <- rbind(df, list(row=r, col=c, bin=b, quality=Q) )
        }
        }
    }

    df$bin <- factor(df$bin)
    levels(df$bin) <- c("Run start", "Run middle", "Run end")
    p5 = ggplot(df, aes(col, row)) + 
        geom_tile(aes(fill = quality), color="white") +
        scale_fill_gradient(low = "white", high="steelblue") +
        xlab("Column") +
        ylab("Row") +
        facet_wrap(~bin, ncol=3) +
        theme_linedraw() +
        theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
        #theme(text = element_text(size = 40), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), legend.text=element_text(size=12))
    ggsave(filename = file.path(output.dir, "flowcell_overview.png"), width = 7, height = 4, dpi=300, plot = p5)
    ggsave(filename = file.path(output.dir, "flowcell_overview.screen.png"), width = 7, height = 4, dpi=90, plot = p5)

    flog.info(paste(sep = "", flowcell, ": plotting cumulative yield summary"))
    p5a = ggplot(d, aes(x=start_time/3600, y=cumulative.bases.time, colour = Q_cutoff)) + 
        geom_line(size = 1) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        xlab("Elapsed time (hrs)") +
        ylab("Total yield in bases") +
        scale_colour_discrete(guide = guide_legend(title = "Reads")) +
        scale_x_continuous() +
        theme_linedraw()
        #theme(text = element_text(size = 15))

    ggsave(filename = file.path(output.dir, "cumulative_yield.png"), width = 7, height = 4, dpi=300, plot = p5a)
    ggsave(filename = file.path(output.dir, "cumulative_yield.screen.png"), width = 7, height = 4, dpi=90, plot = p5a)

    flog.info(paste(sep = "", flowcell, ": plotting flowcell yield summary"))
    p6 = ggplot(d, aes(x=sequence_length_template, y=cumulative.bases, colour = Q_cutoff)) + 
        geom_line(size = 1) + 
        xlab("Minimum read length") +
        ylab("Total yield in bases") +
        scale_colour_discrete(guide = guide_legend(title = "Reads")) +
        #theme(text = element_text(size = 15))
        theme_linedraw()
    xmax = max(d$sequence_length_template[which(d$cumulative.bases > 0.01 * max(d$cumulative.bases))])
    p6 = p6 + scale_x_continuous(limits = c(0, xmax))

    ggsave(filename = file.path(output.dir, "yield_summary.png"), width = 7, height = 4, dpi=300, plot = p6)
    ggsave(filename = file.path(output.dir, "yield_summary.screen.png"), width = 7, height = 4, dpi=90, plot = p6)
    
    flog.info(paste(sep = "", flowcell, ": plotting sequence length over time"))
    e = subset(d, Q_cutoff=="All reads")
    e$Q = paste(">=", q, sep="")
    e$Q[which(e$mean_qscore_template<q)] = paste("<", q, sep="")
    p7 = ggplot(e, aes(x=start_time/3600, y=sequence_length_template, colour = Q, group = Q)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        geom_smooth() + 
        xlab("Elapsed time (hrs)") + 
        ylab("Mean read length") + 
        ylim(0, NA) +
        theme_linedraw()
    ggsave(filename = file.path(output.dir, "length_by_hour.png"), width = 7, height = 4, dpi=300, plot = p7)
    ggsave(filename = file.path(output.dir, "length_by_hour.screen.png"), width = 7, height = 4, dpi=90, plot = p7)
    
    flog.info(paste(sep = "", flowcell, ": plotting Q score over time"))
    p8 = ggplot(e, aes(x=start_time/3600, y=mean_qscore_template, colour = Q, group = Q)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        geom_smooth() + 
        xlab("Elapsed time (hrs)") + 
        ylab("Mean Q score") + 
        ylim(0, NA) +
        theme_linedraw()
    ggsave(filename = file.path(output.dir, "q_by_hour.png"), width = 7, height = 4, dpi=300, plot = p8)
    ggsave(filename = file.path(output.dir, "q_by_hour.screen.png"), width = 7, height = 4, dpi=90, plot = p8)
    
    flog.info(paste(sep = "", flowcell, ": plotting reads per hour"))
    f = d[c("hour", "reads_per_hour", "Q_cutoff")]
    f = f[!duplicated(f),]
    g = subset(f, Q_cutoff=="All reads")
    h = subset(f, Q_cutoff==q_title)
    max = max(f$hour)
    # all of this is just to fill in hours with no reads recorded
    all = 0:max
    add.g = all[which(all %in% g$hour == FALSE)]
    if(length(add.g)>0){
        add.g = data.frame(hour = add.g, reads_per_hour = 0, Q_cutoff = "All reads")
        g = rbind(g, add.g)
    }
    add.h = all[which(all %in% h$hour == FALSE)]
    if(length(add.h)>0){
        add.h = data.frame(hour = add.h, reads_per_hour = 0, Q_cutoff = q_title)
        h = rbind(h, add.h)
    }
    i = rbind(g, h)
    i$Q_cutoff = as.character(i$Q_cutoff)
    i$Q_cutoff[which(i$Q_cutoff==q_title)] = paste("Q>=", q, sep="")
    p9 = ggplot(i, aes(x=hour, y=reads_per_hour, colour = Q_cutoff, group = Q_cutoff)) + 
        geom_vline(xintercept = muxes, colour = 'red', linetype = 'dashed', alpha = 0.5) +
        #geom_point() +
        geom_line(size=1) +
        xlab("Elapsed time (hrs)") + 
        ylab("Number of reads per hour") + 
        ylim(0, NA) + 
        scale_color_discrete(guide = guide_legend(title = "Reads")) +
        theme_linedraw()
    ggsave(filename = file.path(output.dir, "reads_per_hour.png"), width = 7, height = 4, dpi=300, plot = p9)
    ggsave(filename = file.path(output.dir, "reads_per_hour.screen.png"), width = 7, height = 4, dpi=90, plot = p9)
    
    flog.info(paste(sep = "", flowcell, ": plotting read length vs. q score scatterplot"))
    p10 = ggplot(subset(d, Q_cutoff=="All reads"), aes(x = sequence_length_template, y = mean_qscore_template, colour = events_per_base)) + 
        geom_point(alpha=0.05, size = 0.4) + 
        scale_x_log10(minor_breaks=log10_minor_break()) + 
        labs(colour='Events per base\n(log scale)\n')  + 
        theme(text = element_text(size = 15)) +
        xlab("Read length") +
        ylab("Mean Q score of read") +
        theme_linedraw()

    if(max(d$events_per_base, na.rm=T)>0){
        # a catch for 1D2 runs which don't have events per base
        p10 = p10 + scale_colour_distiller(trans = "log", labels = scientific, palette = 'Spectral') 
    }

    ggsave(filename = file.path(output.dir, "length_vs_q.png"), width = 7, height = 4, dpi=300, plot = p10)
    ggsave(filename = file.path(output.dir, "length_vs_q.screen.png"), width = 7, height = 4, dpi=90, plot = p10)
    
    return(d)
}

if(file_test("-f", input.file)==TRUE){
    # if it's an existing file (not a folder) just run one analysis
    d = single.flowcell(input.file, output.dir, q)
}else{
    #WTF
    flog.warn(paste("Could find a sequencing summary file in your input which was: ", 
        input.file, 
        "\nThe input must be a sequencing_summary.txt file produced by Albacore"))
}
