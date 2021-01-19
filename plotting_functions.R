readgwas <- function(startpath = "/zone_gwas/output", suffix = "*.assoc.txt"){
  gwas = list.files(path = paste(here::here(), startpath, sep = ""),
             pattern = glob2rx(suffix), full.names = TRUE) %>%
    purrr::map(~read_tsv(.x, col_names = TRUE, col_types=cols(rs = col_character()))) %>%
    purrr::reduce(bind_rows) %>%
    mutate(q = qvalue(p_score)$qvalues) %>% 
    mutate(chr = as.numeric(str_split_fixed(rs, ":", n = 2)[,1])) %>%
    mutate(pos = as.numeric(str_split_fixed(rs, ":", n = 2)[,2])) %>%
    #left_join(., genes(.)) %>%
    #rename(p = p_score) %>% 
    select(rs, chr, pos ,af, p = p_score, q)
    
    don = gwas %>% 
      group_by(chr) %>% 
      summarise(chr_len=max(pos)) %>% 
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%
      left_join(gwas, ., by=c("chr"="chr")) %>%
      arrange(chr, pos) %>%
      mutate( BPcum=pos+tot) %>% 
      ungroup()
  return(don)
}

readgwas2 <- function(filepath){
    gwas = 
      read_tsv(filepath, col_names = TRUE, col_types=cols(rs = col_character())) %>% 
    mutate(q = qvalue(p_score)$qvalues) %>% mutate(chr = as.numeric(str_split_fixed(rs, ":", n = 2)[,1])) %>% 
    mutate(pos = as.numeric(str_split_fixed(rs, ":", n = 2)[,2])) %>% 
    #left_join(., genefromsnp(.)) %>% 
    select(rs, chr, pos, beta , af, p = p_score, q)
  
  don = gwas %>% 
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by=c("chr"="chr")) %>%
    arrange(chr, pos) %>%
    mutate( BPcum=pos+tot) %>% 
    ungroup()
  
  return(don)
}


ggmanhattan = function(inputfile, 
                       prune = 1.0, 
                       value = p, 
                       alpha = 0.5, 
                       pcol = "p_score", 
                       pointsize = 1.0,
                       colors = c("grey10", "grey55"), 
                       sigsnps = NULL, 
                       sigsnps2 = NULL,
                       sigsnpcolor = "springgreen3",
                       sigsnpcolor2 = "goldenrod1"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(value)
  gwas = inputfile %>% #reads in input file
    #select(rs, p = pcol) %>% 
    #mutate(chr = case_when(is.null(chr) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,1]),
    #                       TRUE ~ chr)) %>% 
    #mutate(pos = case_when(is.null(pos) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,2]),
    #                       TRUE ~ pos)) %>% 
    #mutate(chr = as.numeric(str_split_fixed(rs, ":", n = 2)[,1]))%>% #extracts information from SNP name which should be chr:pos
    #mutate(pos = as.numeric(str_split_fixed(rs, ":", n = 2)[,2])) %>% #This holds as long as this was calculated in my imp pipeline
    #mutate(q = qvalue(p)$qvalues) %>% #transforms p-values to q-values
    #ifelse(v == "p", filter(p<prune), filter(q<prune))
    filter(., if (v == "p") p < prune else q < prune)

  
  axisdf = gwas %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )  
  
  #newsigsnps = left_join(sigsnps, don, by = "rs")
  #print(head(newsigsnps))
  #v from above 
  gwas_plot = ggplot(gwas, aes(x=BPcum, y=-log10(!!v))) + 
    #geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    geom_point(aes(color=as.factor(chr)), alpha=alpha, size = pointsize) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$chr[c(TRUE, FALSE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
          y = case_when(v == "p" ~ expression(paste(-log10, "(", italic('p'), ")")),
                        v == "q" ~ expression(paste(-log10, "(", italic('q'), ")"))))+
    theme_bw() +
    # geom_point(data = subset(gwas, rs %in% sigsnps), color=sigsnpcolor, size=1.3) +
    # #geom_point(data = newsigsnps, aes(x=BPcum, y=-log10(!!v), color=sigsnpcolor), size=1.3) +
    # geom_point(data = subset(gwas, rs %in% sigsnps2), color=sigsnpcolor2, size=1.3) +
    geom_point(data = subset(gwas, rs %in% sigsnps), color=sigsnpcolor, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(gwas, rs %in% sigsnps2), color=sigsnpcolor2, size = pointsize, alpha = 0.5) +
    geom_hline(yintercept = case_when(v == "p" ~ -log10(1e-5),
                                      v == "q" ~ -log10(0.1)), 
               color = "red", 
               size = 0.25) +
    theme(legend.position="none", 
          panel.border = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )#+
    # theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
    #       axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),  
    #       axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
    #       axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}




ggmanhattan2 = function(inputfile, 
                       prune = 1.0, 
                       value = p, 
                       alpha = 0.5, 
                       pcol = "p_score", 
                       pointsize = 1.0,
                       colors = c("grey10", "grey55"), 
                       sigsnps = NULL, 
                       sigsnps2 = NULL,
                       sigsnps3 = NULL,
                       sigsnpcolor = "springgreen3",
                       sigsnpcolor2 = "goldenrod1",
                       sigsnpcolor3 = "royalblue4"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(value)
  #print(v)
  gwas = inputfile %>% #reads in input file
    select(rs, p = pcol) %>% 
    # mutate(chr = case_when(is.null(chr) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,1]),
    #                        TRUE ~ chr)) %>%
    # mutate(pos = case_when(is.null(pos) ~ as.integer(str_split_fixed(snp, ":", n = 2)[,2]),
    #                        TRUE ~ pos)) %>%
    mutate(chr = as.numeric(str_split_fixed(rs, ":", n = 2)[,1]))%>% #extracts information from SNP name which should be chr:pos
    mutate(pos = as.numeric(str_split_fixed(rs, ":", n = 2)[,2])) %>% #This holds as long as this was calculated in my imp pipeline
    mutate(q = qvalue(p)$qvalues) %>% #transforms p-values to q-values
    #ifelse(v == "p", filter(p<prune), filter(q<prune))
    filter(., if (v == "p") p < prune else q < prune)
  
  don = gwas %>% 
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by=c("chr"="chr")) %>%
    arrange(chr, pos) %>%
    mutate(BPcum=pos+tot) %>% 
    ungroup()
  
  axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  #newsigsnps = left_join(sigsnps, don, by = "rs")
  #print(head(newsigsnps))
  #v from above 
  gwas_plot = 
    don %>% 
    filter(!rs %in% c(sigsnps, sigsnps2, sigsnps3)) %>% 
    ggplot(
    aes(x=BPcum, y=-log10(!!v))) + 
    #geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    geom_point(aes(color=as.factor(chr)), alpha=alpha, size = pointsize) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$chr[c(TRUE, FALSE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
         y = case_when(v == "p" ~ expression(paste(-log10, "(", italic('p'), ")")),
                       v == "q" ~ expression(paste(-log10, "(", italic('q'), ")"))))+
    theme_bw() +
    # geom_point(data = subset(gwas, rs %in% sigsnps), color=sigsnpcolor, size=1.3) +
    # #geom_point(data = newsigsnps, aes(x=BPcum, y=-log10(!!v), color=sigsnpcolor), size=1.3) +
    # geom_point(data = subset(gwas, rs %in% sigsnps2), color=sigsnpcolor2, size=1.3) +
    geom_point(data = subset(don, rs %in% sigsnps), color=sigsnpcolor, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, rs %in% sigsnps2), color=sigsnpcolor2, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, rs %in% sigsnps3), color=sigsnpcolor3, size = pointsize, alpha = 0.5) +
    geom_hline(yintercept = case_when(v == "p" ~ -log10(1e-5),
                                      v == "q" ~ -log10(0.1)), 
               color = "red", 
               size = 0.25) +
    theme(legend.position="none", 
          panel.border = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )#+
  # theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),  
  #       axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
  #       axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}


selscan_manhattans = function(inputfile,
                              col = ihs,
                              alpha = 0.5,
                              pointsize = 1.0,
                              colors = c("grey10", "grey55"),
                              ylab = "",
                              sigsnps = NULL,
                              sigsnps2 = NULL,
                              sigsnps3 = NULL,
                              sigsnpcolor = "springgreen3",
                              sigsnpcolor2 = "goldenrod1",
                              sigsnpcolor3 = "royalblue4"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(col)
  #print(v)
  gwas = inputfile
  
  don = gwas %>%
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    ungroup()
  
  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  gwas_plot =
    don %>%
    ggplot(
      aes(x=BPcum, y=!!v)) +
    #geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    geom_point(aes(color=as.factor(CHR)), alpha=alpha, size = pointsize) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$CHR[c(TRUE, FALSE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
         y = "")+
    theme_bw() +
    geom_point(data = subset(don, SNP %in% sigsnps), color=sigsnpcolor, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, SNP %in% sigsnps2), color=sigsnpcolor2, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, SNP %in% sigsnps3), color=sigsnpcolor3, size = pointsize, alpha = 0.5) +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )#+
  # theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),
  #       axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
  #       axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}



readevec <- function(evecfile, program = "smartpca"){
  if (program == "smartpca"){ 
    evec = read_delim(evecfile, 
                    col_names = c("identifier", paste("pc", seq(1,10), sep = ""), "pop"), 
                    delim = " ", trim_ws = TRUE) %>% 
      separate("identifier", c("family", "id"),sep = ":")
    }
  else{
    evec = read_delim(evecfile, 
                      col_names = c("family", "identifier", paste("pc", seq(1,20), sep = "")), 
                      delim = " ", trim_ws = TRUE)
  }
  return(evec)
}

ggpca <- function(evec, pc_a=pc1, pc_b=pc2, method = "smartpca", size = 2, evalfile = NULL, colorby = NULL, title =NULL, legend = FALSE, alpha = 0.5){
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  require(readr)
  a = enexpr(pc_a)
  b = enexpr(pc_b)
  color = enexpr(colorby)
  if (method == "smartpca"){
  eval = slice(evec, 1) %>% 
    select(pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)
  evec = slice(evec, 2:n()) 
  eval_a = round(eval[as.integer(str_replace(a, "pc", "")),]$eval/sum(eval)*100, digits = 2)
  eval_b = round(eval[1,as.integer(str_replace(b, "pc", ""))]/sum(eval)*100, digits = 2)
  }
  if (method == "plink"){
    eval = read_table(evalfile, col_names = "eval")
    eval_a = round(eval[as.integer(str_replace(a, "pc", "")),]$eval/sum(eval)*100, digits = 2)
    eval_b = round(eval[as.integer(str_replace(b, "pc", "")),]$eval/sum(eval)*100, digits = 2)
  }
  colour = c("3"="springgreen3", "9"="slateblue2", "1"="tomato2", "5"="goldenrod1", "6"="gray50", "8"="gray17", "4"="brown", "2"="darkslategray4", "7"="deeppink3")
  region_names = c("1" = "Desert", "2" = "Southeast", "3" = "High Plains", "4" = "Rainforest", "5" = "Arid Prairie", "6" = "Foothills", "7" = "Forested\nMountains", "8" = "Fescue Belt", "9" = "Upper Midwest &\nNortheast")
  pca_plot = ggplot(evec, aes(x=!!a, y=!!b, color = as.factor(!!color)))+
    geom_point(alpha = alpha, show.legend = legend, size = size)+
    ggtitle(if (is.null(title)) "" else paste(title," (n = ", count(evec)$n, " individuals)", sep = ""))+ 
    xlab(paste("Principal Component ", str_replace(a, "pc", ""), "\n", eval_a, " % of genetic variance explained"))+
    ylab(paste("Principal Component ", str_replace(b, "pc", ""),", \n", eval_b, " % of genetic variance explained"))+
    scale_color_manual("Ecoregion", values=colour, labels = region_names[evec$zone])
  return(pca_plot)
}

mapplotting <- function(data, colorby = "zone", title = NULL, legend_title = "", vcolor = NULL, max = 4, min = 0.1, alpha = 0.6){
  colour = c("3"="springgreen3", "9"="slateblue2", "1"="tomato2", "5"="goldenrod1", "6"="gray50", "8"="gray17", "4"="brown", "2"="darkslategray4", "7"="deeppink3")
  all_states<-map_data("state")
  region_names = c("1" = "Desert", "2" = "Southeast", "3" = "High Plains", "4" = "Rainforest", "5" = "Arid Prairie", "6" = "Foothills", "7" = "Forested\nMountains", "8" = "Fescue Belt", "9" = "Upper Midwest &\nNortheast")
  #colorby <- enquo(colorby)
  if (colorby == "zone"){
  mapplot = ggplot()+
    geom_polygon(data=all_states, aes(x=long, y=lat, group=group), color="black", fill="white")+
    geom_count(aes(x = data$x, y = data$y, colour=factor(data$zone)), alpha = .6)+
    scale_size(range = c(min, max), guide = FALSE)+ #This changes max dot size for altering figure size
    scale_color_manual("Climate Zone", values=colour, labels = paste("Zone ", levels(data$zone), " : ", table(data$zone), sep = ""))+
    ggtitle(title)+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme_map()
  }
  if (colorby != "zone"){
  #colorby <- enquo(colorby)
    if (colorby == "meantemp"){vcolor = "B" 
                              dir = 1}
    if (colorby == "precip"){vcolor = "D" 
                            dir = -1}
    if (colorby == "elev"){vcolor = "B" 
                            dir = 1}
  mapplot = ggplot()+
    geom_polygon(data=all_states, aes(x=long, y=lat, group=group), color="black", fill="white")+
    geom_count(data = data, aes(x = x, y = y, color=!!sym(colorby)), alpha = 0.6)+
    scale_size(range = c(min, max), guide = FALSE)+ #This changes max dot size for altering figure size
    scale_color_viridis(option = vcolor, direction = -1)+
    labs(color = legend_title)+
    ggtitle(title)+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme_map()
  }
return(mapplot)
}




mapplotting2 <- function(data, title = "", legend_title = "", vcolor = NULL, max = 4, min = 0.1, alpha = 0.6, facet = NULL){
  colour = c("3"="springgreen3", "9"="slateblue2", "1"="tomato2", "5"="goldenrod1", "6"="gray50", "8"="gray17", "4"="brown", "2"="darkslategray4", "7"="deeppink3")
  region_names = c("1" = "Desert", "2" = "Southeast", "3" = "High Plains", "4" = "Rainforest", "5" = "Arid Prairie", "6" = "Foothills", "7" = "Forested\nMountains", "8" = "Fescue Belt", "9" = "Upper Midwest &\nNortheast")
  all_states<-map_data("state")
    mapplot = ggplot()+
      geom_count(aes(x = data$x, y = data$y, colour=factor(data$zone)), alpha = alpha)+
      scale_size(range = c(min, max), guide = FALSE)+ #This changes max dot size for altering figure size
      scale_color_manual("Climate Zone", values=colour, labels = region_names[data$zone])+
      #scale_color_manual("Climate Zone", values=colour, labels = paste("Zone ", levels(data$zone), " : ", table(data$zone), sep = ""))+
      ggtitle(title)+
      geom_polygon(data = map_data("state"), aes(x = long, y = lat, group = group), color = "black", fill = NA)+
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      theme_map()
  return(mapplot)
}

locustimeplot <- function(stratfile, n = 1, linesize = 5, linecolor = "firebrickred", title = NULL, error = FALSE, flip = "no"){
  nrows = n * 24
  aftime = read_table(stratfile, col_names = TRUE, col_types = cols(SNP = col_character())) %>% 
  mutate(CLST = case_when(CLST == "75to85" ~ 1980,
                          CLST == "86to95" ~ 1990,
                          TRUE ~ as.numeric(CLST))) %>% 
  mutate(AF = MAC/NCHROBS) %>%
  mutate(AF = case_when(flip == "yes" ~ 1-AF,
                        TRUE ~ AF)) %>% 
    left_join(pb_gpsm, by = c("SNP" = "rs")) %>% 
    top_n(nrows, wt = -log10(p_score))
  locusplot =
    ggplot(aftime, aes(CLST, AF, color = SNP))+
    geom_smooth(size = linesize, se = error, color = linecolor)+
    labs(x = "Year of birth", y = "AF")+
    ggtitle(title)+
    theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "bold"))
return(locusplot)
}

ggqq <- function(pvector){
  pvector = pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  pdf = data.frame(observed = -log10(sort(pvector,decreasing=FALSE)), expected = -log10(ppoints(length(pvector))))
  qqplotted = ggplot(pdf, aes(expected, observed))+
    geom_point()+
    geom_abline(intercept = 0,
                slope = 1,
                colour = "red")+
    labs(x = expression(paste("Expected ", -log10, "(", italic('p'), ")")),
         y = expression(paste("Observed ", -log10, "(", italic('p'), ")")))
  return(qqplotted)
}  

ggqq <- function(pvector, size = NULL){
  pvector = pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  pdf = data.frame(observed = -log10(sort(pvector,decreasing=FALSE)), expected = -log10(ppoints(length(pvector))))
  qqplotted = ggplot(pdf, aes(expected, observed))+
    geom_point(size = size)+
    geom_abline(intercept = 0,
                slope = 1,
                colour = "red")+
    labs(x = expression(paste("Expected ", -log10, "(", italic('p'), ")")),
         y = expression(paste("Observed ", -log10, "(", italic('p'), ")")))+
    theme_cowplot()
  return(qqplotted)
}

mytheme = 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 33, hjust = 1),
        strip.text.x = element_text(size = 8),
        legend.text = element_text(size = 10),
        legend.title = element_blank())


