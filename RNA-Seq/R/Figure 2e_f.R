
##A example of genes that their expression increased or decreased monotonously at different stages of muscle development
library(dplyr)
library(stringr)
library(rlang)

OriExp <- read.csv("F:\\Manuscript\\RNA\\Ori\\data\\PCA\\featureCount_result\\GeneExpressionAllSamples.tpm.csv", header = T, row.names = 1)

OriExp_f <- select(OriExp, !contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1")))

OriExpMu <- select(OriExp_f, contains("Mu"))

##ENSSSCG00000023128 (Gene expression increased monotonously) IGF2BP1
##ENSSSCG00000038323 (Gene expression decreased monotonously) ACTN3

##IGF2BP1
gene <- c("ENSSSCG00000023128")

OriExpMugene <- filter(OriExpMu, rownames(OriExpMu) %in% gene)

F40 <- OriExpMugene %>% select(.,contains("40")) %>% t() %>% as.data.frame()
F70 <- OriExpMugene %>% select(.,contains("70")) %>% t() %>% as.data.frame()
D1 <- OriExpMugene %>% select(.,contains("115")) %>% t() %>% as.data.frame()
D168 <- OriExpMugene %>% select(.,contains("168")) %>% t() %>% as.data.frame()

StgCombine <- rbind(F40,F70,D1,D168)

tmp <- unlist(str_split(rownames(StgCombine),"_"))

stage <- str_sub(tmp[seq(1,length(tmp),3)],3)
stage <- replace(stage,grep("40",stage),"F40") %>% replace(.,grep("70",stage),"F70") %>% replace(.,grep("115",stage),"D1") %>% replace(.,grep("168",stage),"D168")

stage <- factor(stage, levels= c("F40", "F70", "D1", "D168"))

StgCombine$Stage <- stage

datPlot <- StgCombine
colnames(datPlot) <- c("Gene","Stage")

library(ggplot2)
library(ggsci)

p <- ggplot(datPlot, aes(x=Stage,y=log2(Gene+1),fill=Stage)) + geom_bar(stat = "summary", fun=mean, width = 0.5, color="black") + stat_summary(fun.data = 'mean_sdl', geom="uperrorbar", color="black", width=0.2, position = position_dodge(.9))+theme_bw() + scale_fill_lancet() + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("IGF2BP1 expression log2(TPM)")+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) +theme(panel.border = element_blank(),axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"))

pdf("IGF2BP1.diffStg.Exp.pdf", width = 4, height = 4)
p
dev.off()


##ACTN3
gene <- c("ENSSSCG00000038323")

OriExpMugene <- filter(OriExpMu, rownames(OriExpMu) %in% gene)

F40 <- OriExpMugene %>% select(.,contains("40")) %>% t() %>% as.data.frame()
F70 <- OriExpMugene %>% select(.,contains("70")) %>% t() %>% as.data.frame()
D1 <- OriExpMugene %>% select(.,contains("115")) %>% t() %>% as.data.frame()
D168 <- OriExpMugene %>% select(.,contains("168")) %>% t() %>% as.data.frame()

StgCombine <- rbind(F40,F70,D1,D168)

tmp <- unlist(str_split(rownames(StgCombine),"_"))

stage <- str_sub(tmp[seq(1,length(tmp),3)],3)
stage <- replace(stage,grep("40",stage),"F40") %>% replace(.,grep("70",stage),"F70") %>% replace(.,grep("115",stage),"D1") %>% replace(.,grep("168",stage),"D168")

stage <- factor(stage, levels= c("F40", "F70", "D1", "D168"))

StgCombine$Stage <- stage

datPlot <- StgCombine
colnames(datPlot) <- c("Gene","Stage")

library(ggplot2)
library(ggsci)

p <- ggplot(datPlot, aes(x=Stage,y=log2(Gene+1),fill=Stage)) + geom_bar(stat = "summary", fun=mean, width = 0.5, color="black") + stat_summary(fun.data = 'mean_sdl', geom="uperrorbar", color="black", width=0.2, position = position_dodge(.9))+theme_bw() + scale_fill_lancet() + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("ACTN3 expression log2(TPM)")+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16)) +theme(panel.border = element_blank(),axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"))

pdf("ACTN3.diffStg.Exp.pdf", width = 4, height = 4)
p
dev.off()


#---------------------------------------------------
##A R function to keep the unilateral error line
geom_uperrorbar <- function(mapping = NULL, data = NULL,
                          stat = "identity", position = "identity",
                          ...,
                          na.rm = FALSE,
                          orientation = NA,
                          show.legend = NA,
                          inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
  default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
    alpha = NA),

  draw_key = draw_key_path,

  required_aes = c("x|y", "ymin|xmin", "ymax|xmax"),

  setup_params = function(data, params) {
    GeomLinerange$setup_params(data, params)
  },

  extra_params = c("na.rm", "orientation"),

  setup_data = function(data, params) {
    data$flipped_aes <- params$flipped_aes
    data <- flip_data(data, params$flipped_aes)
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)
    data <- transform(data,
      xmin = x - width / 2, xmax = x + width / 2, width = NULL
    )
    flip_data(data, params$flipped_aes)
  },

  draw_panel = function(data, panel_params, coord, width = NULL, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    #x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x,    NA, data$xmin, data$xmax))
    #y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$ymin, NA, data$ymin, data$ymin))
    sel <- data$y < 0 
    data$ymax[sel] <- data$ymin[sel]
    x <- as.vector(rbind(data$xmin, data$xmax, NA, data$x,    data$x))
    y <- as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y))
    data <- new_data_frame(list(
      x = x,
      y = y,
      colour = rep(data$colour, each = 5),
      alpha = rep(data$alpha, each = 5),
      size = rep(data$size, each = 5),
      linetype = rep(data$linetype, each = 5),
      group = rep(1:(nrow(data)), each = 5),
      row.names = 1:(nrow(data) * 5)
    ))
    data <- flip_data(data, flipped_aes)
    GeomPath$draw_panel(data, panel_params, coord)
  }
)

new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }
  
  class(x) <- "data.frame"
  
  attr(x, "row.names") <- .set_row_names(n)
  x
}