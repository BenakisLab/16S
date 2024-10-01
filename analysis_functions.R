if(!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(
  BiocManager,
  tidyverse,
  ggpubr,
  ggsci,
  rstatix,
  picante,
  phyloseq,
  vegan,
  ANCOMBC,
  phangorn,
  GUniFrac,
  zoo
)

######################### DATA WRANGLING #########################

taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

meta_to_df <- function(ps) {
  return(as(sample_data(ps), "data.frame"))
}

ps_to_asvtab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}


#########################    IMPORTING DATA #########################
split_taxonomy <- function(asvtab) {
  # select taxa column and split taxonomy column into separate columns using ; delimiter
  taxonomy_phylo <- dplyr::select(asvtab, "taxonomy") %>%
    separate("taxonomy",
             c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";")
  return(taxonomy_phylo)
}

format_taxonomy <- function(ps) {
  ps_tax <- taxonomy(ps) %>%
    mutate_all(na_if, "")
  
  ps_tax[] <- t(na.locf(t(ps_tax))) %>%
    as.data.frame()
  # add unknown to genera which match family
  ps_tax <- ps_tax %>%
    mutate(Genus = case_when(Genus == Family ~ paste0("unknown_", Genus),
                             TRUE ~ Genus)) %>%
    mutate(Family = case_when(Family == Order ~ paste0("unknown_", Family),
                              TRUE ~ Family)) %>%
    mutate(Order = case_when(Order == Class ~ paste0("unknown_", Order),
                             TRUE ~ Order)) %>%
    mutate(Class = case_when(Class == Phylum ~ paste0("unknown_", Class),
                             TRUE ~ Class)) %>%
    mutate(Phylum = case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                              TRUE ~ Phylum)) %>%
    mutate(ASV = rownames(ps_tax)) %>%
    mutate(Highest_classified = paste(rownames(ps_tax), Genus, sep = "; "))
  
  rownames(ps_tax) <- ps_tax[, "Highest_classified"]
  taxa_names(ps) <- ps_tax[, "Highest_classified"]
  tax_table(ps) <- tax_table(as.matrix(ps_tax))
  
  return(ps)
}


read_tab_delim <- function(df) {
  # read all tab delimited files using these params
  df_out <-
    read.table(
      df,
      row.names = 1,
      header = 1,
      sep = "\t",
      comment.char = "",
      check.names = F
    )
  return(df_out)
}

load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  # convert to phyloseq and return list
  phylo_asv <- otu_table(asvtab, taxa_are_rows = T)
  
  phylo_tax <- tax_table(as.matrix(taxa))
  
  phylo_map <- sample_data(mapping)
  
  if (exists("tree")) {
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}

import_as_pseq <- function(asvtab, mapping, tree = NULL) {
  # read files
  asvtab_taxa <- read_tab_delim(asvtab)
  metadata <- read_tab_delim(mapping)
  # convert otu into matrix and drop taxonomy col
  asv_matrix <- as.matrix(subset(asvtab_taxa, select = -taxonomy))
  # split taxonomy
  taxonomy_phylo <- split_taxonomy(asvtab_taxa)
  # make sure taxa table row names match otu
  rownames(taxonomy_phylo) <- rownames(asvtab_taxa)
  # use load phylo function to convert all files to phyloseq objs
  if (exists("tree")) {
    # if user wants to input a phylogenetic tree this code is called
    out <-
      load_phylo(asv_matrix, taxonomy_phylo, metadata, tree = tree)
  }
  else {
    # without a tree
    out <- load_phylo(asv_matrix, taxonomy_phylo, metadata)
  }
  return(out)
}

######################### NORMALISATION #########################

# normalise data with minimum sum scaling
transform <- function(ps, transform = "mss") {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_asvtab(ps)
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }
  
  if (transform %in% c("mss", "relative")) {
    if (transform == "mss") {
      ps_t <- t(min(colSums(x)) * t(x) / colSums(x))
      
    } else if (transform == "relative") {
      ps_t <- t(100 * t(x) / colSums(x))
      
    }
    otu_table(ps)@.Data <- ps_t
    
    return(ps)
    
  } else {
    print("Not a valid transform, exiting")
    stop()
  }
  
}


######################### ALPHA DIVERSITY #########################

# calculate shannon effective
Shannon.E <- function(x) {
  shannon.e <- apply(x, 1, function(x) round(exp(-sum(x[x > 0] / sum(x) * log(x[x > 0] / sum(x)))), digits = 2))
  return(shannon.e)
}

Simpson <- function(x) {
  simpson <- vegan::diversity(x, index = "invsimpson")
  return(simpson)
}

# calculate faiths phylogenetic diversity
Faiths <- function(x, tree) {
  pd <- pd(x, tree, include.root = F)
  return(pd[, 1])
}

# calculate richness
Richness <- function(x, detection = 0.5) {
  observed <- apply(x > detection, 1, sum)
  return(observed)
}

# calculate all alpha diversity matrices and return dataframe
calc_alpha <- function(ps, metrics=c("Richness", "Shannon.Effective", "Faiths.PD", "Inverse.Simpson"), ...) {
  mat_in <- ps_to_asvtab(ps) %>% 
    t()
  
  tree <- phy_tree(ps)
  
  diversity <-
    setNames(data.frame(matrix(ncol = length(metrics), nrow = nsamples(ps))),
             metrics)
  rownames(diversity) <- rownames(meta_to_df(ps)) 
  
  if ("Richness" %in% metrics){
    diversity$Richness <- Richness(mat_in)
  }
  if ("Shannon.Effective" %in% metrics){
    diversity$Shannon.Effective <- Shannon.E(mat_in)
  }
  if ("Faiths.PD" %in% metrics){
    diversity$Faiths.PD <- Faiths(mat_in, tree)
  }
  if ("Inverse.Simpson" %in% metrics){
    diversity$Inverse.Simpson <- Simpson(mat_in)
  }
  
  # merge alpha diversity df and metadata by rownames (0)
  diversity <- merge(diversity, meta_to_df(ps), by=0)

  return(diversity)
}

######################### BETA DIVERSITY #########################

# Calculate Gunifrac in phyloseq objects
phyloseq_gunifrac <- function(ps, asdist = TRUE) {
  # input phyloseq obj
  # returns dist matrix of GUnifrac distance which can
  # be used as input for ordinate() function
  
  #extract otu table and tree
  asvtab <- t(ps_to_asvtab(ps))
  tree <- phy_tree(ps)
  
  # root tree at midpoint
  rooted_tree <- midpoint(tree)
  
  # calculate gunifrac
  gunifrac <- GUniFrac(asvtab, rooted_tree,
                       alpha = c(0.0, 0.5, 1.0))$unifracs
  # alpha param is weight on abundant lineages, 0.5 has best power
  # so we will extract this
  if (asdist == T) {
    return(stats::as.dist(gunifrac[, , "d_0.5"]))
  }
  else{
    return(gunifrac[, , "d_0.5"])
  }
}


calc_betadiv <- function(ps, dist, ord_method = "NMDS") {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {
    
    if (dist %in% c("unifrac", "wunifrac", "bray")) {
      dist_metric <- distance(ps, dist)
      ord <- ordinate(ps, ord_method, dist_metric)
      return_list <- list("Distance_Matrix" = dist_metric,
                          "Ordination" = ord)
      return(return_list)
    }
    else if (dist %in% c("gunifrac")) {
      gu <- phyloseq_gunifrac(ps)
      ord <- ordinate(ps, ord_method, gu)
      
      return_list <- list("Distance_Matrix" = gu,
                          "Ordination" = ord)
      return(return_list)
    }
    else {
      print(
        "Distance metric not supported, supported metrics are; bray, unifrac, wunifrac, gunifrac"
      )
    }
  }
  else {
    print("Ordination method not supported, supported methods are: NMDS, MDS, PCoA")
    stop()
  }
  
}


######################### Differential Abundance #########################
add_taxonomy_da <- function(ps, list_da_asvs, da_res) {
  taxonomy <- as.data.frame(tax_table(ps))
  
  da_taxonomy <- taxonomy %>%
    select(Highest_classified) %>%
    filter(Highest_classified %in% list_da_asvs)
  
  da_tax_out <- da_res %>%
    rownames_to_column(var = "Group") %>%
    column_to_rownames("ASV") %>%
    merge(da_taxonomy, by = 0)
  
  return(da_tax_out)
}

filt_low_variance <- function(ps, var, thresh=0) {
  # calculate per group abundance and return those higher than threshold 
  # this currently does nothing
  mean_abund <- psmelt(ps) %>% 
    group_by_at(vars(OTU, var)) %>% 
    summarise(mean = mean(Abundance)) %>% 
    filter(all(mean > thresh))
  
  return(mean_abund$OTU)
}

ancom_da <- function(ps,
                     formula,
                     group,
                     ord = NULL,
                     zero_thresh = 0.33,
                     level = NULL,
                     filt_low_var=F,
                     mean_abund_thresh=0,
                     ...) {
  if (!is.null(ord)) {
    sample_data(ps)[[group]] <-
      factor(sample_data(ps)[[group]], levels = ord)
  } else {
    ord <- sort(unique(sample_data(ps)[[group]])) 
  }
  
  levels <-
    c("Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species")
  
  # if level specified check it is valid
  if (!is.null(level)) {
    if (!level %in% levels) {
      stop(
        paste0(
          "Not a valid taxonomic rank, please specify from:",
          levels,
          ", or leave blank"
        )
      )
    }
  }
  
  ps_f <- format_taxonomy(ps)
  
  if (filt_low_var == T){
    low_var_rem <- filt_low_variance(ps_f, group, thresh = 0)
    ps_f <- prune_taxa(low_var_rem, x = ps_f)
  }
  
  res <-
    ancombc2(
      data = ps_f,
      tax_level = level,
      fix_formula = formula,
      p_adj_method = "BH",
      prv_cut = zero_thresh,
      group = group,
      struc_zero = TRUE,
      neg_lb = FALSE,
      global = FALSE,
      ...
    )
  
  # add support for multiple groups here with global and pairwise tests
  res_df <- data.frame(
    # temp fix for ancom version differences
    ASV = unlist(res$res$taxon),
    lfc = unlist(res$res[3]),
    se = unlist(res$res[5]),
    W = unlist(res$res[7]),
    pval = unlist(res$res[9]),
    qval = unlist(res$res[11]),
    da = unlist(res$res[13])
  ) %>%
    mutate(Log2FC = log2(exp(lfc)))
  
  res_da <- res_df %>%
    filter(da == T)
  
  da_asvs <- res_da$ASV
  da_tax <- add_taxonomy_da(ps_f, da_asvs, res_da)
  
  # add group order for easier interpretation
  reference <- paste(ord[1], "vs", ord[2], sep="_")
  da_tax <- da_tax %>% 
    mutate(Group_ord = reference)
  
  return(da_tax)
  
  
}

######################### PLOTTING AND STATS #########################

# calculate colour palettes
calc_pal <- function(ps, group_variable) {
  # update function to accept nonps objects
  meta <- meta_to_df(ps)
  groups <- unique(meta[, group_variable])
  
  if (length(groups) < 10) {
    pal <- pal_npg()(length(groups))
  }
  else {
    print("Exceeded colour palette limit, use custom palette")
  }
}

# calculate adonis R2 and p-value
phyloseq_adonis <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    # convert distance matix object to string
    dist_str <- deparse(substitute(dist_matrix))
    # define formula
    form <- as.formula(paste(dist_str, group_variable, sep="~"))
    ps_ad <- adonis2(form, data = meta_df, ...)
    return(ps_ad)
  }
}


# create statistical formula
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

# check sample sizes for paired test
check_equal_sample_sizes <- function(df, variable_col) {
  # Count the number of occurrences for each group
  sample_sizes <- df %>%
    group_by_at(variable_col) %>%
    summarise(count = n()) %>%
    pull(count)
  
  # Check if all counts are equal
  if (length(unique(sample_sizes)) == 1) {
    return(TRUE)  # All groups have equal sample sizes
  } else {
    return(FALSE) # Sample sizes are not equal across groups
  }
}

# 
stat_plot <- function(df, formula, variable_col, assume_normality=F, multiple_groups=F, paired=F) {
  
  # Check paired
  if (paired == TRUE){
    if (check_equal_sample_sizes(df, variable_col) == F){
      stop("Sample sizes are not equal for paired statistical test")
    }
  }
  
  if (multiple_groups == TRUE) {
    if (paired == TRUE) {
      
      stat_variance <- df %>%
        friedman_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
    else {
      stat_variance <- df %>%
        kruskal_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(formula,
                             comparisons = comparisons_list,
                             p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }
    else {
      stat_test <- df %>%
        wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }
    
  }
}

# plot boxplot with stats
plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group.order = NULL,
                         paired = FALSE,
                         normal=FALSE,
                         ...) {
  # extend color palette with transparent value - required due to way we are
  # layering plot
  if (is.null(cols)) {
    cols <- pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")
  
  if (!is.null(group.order)) {
    df[, variable_col] <-
      factor(df[, variable_col], levels = group.order)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == T) {
    stat_variance <- stat_plot(
      df,
      formula,
      variable_col,
      assume_normality = normal,
      multiple_groups = multiple_groups,
      paired = paired
    )
  }
  stat_test <- stat_plot(
    df,
    formula,
    variable_col,
    assume_normality = normal,
    multiple_groups = multiple_groups,
    paired = paired
  )
  
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(
    df,
    aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    geom_point(size = 1.5, position = position_jitterdodge()) +
    labs(x = xlab, y = ylab) +
    stat_boxplot(color = "black",
                 geom = "errorbar",
                 width = 0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot +
    theme_classic2() +
    ggtitle(p_title) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    rotate_x_text(angle = 45) +
    lemon::coord_capped_cart()
  
  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  }
  else {
    if (multiple_groups == T) {
      label <- "p.adj.signif"
    } else {
      label <- "p.signif"
    }
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = label,
          hide.ns = T,
          inherit.aes = FALSE,
          ...
        )
    }
  
  return(plot_out)
}

# plot scatter plot with correlation if desired
plot_scatter <- function(df,
                         x,
                         y,
                         point_color,
                         line_color,
                         fill_color,
                         xlab,
                         ylab,
                         corr.method = NULL,
                         ...) {
  p <-
    ggplot(data = df, mapping = aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(color = point_color), size = 2.5) +
    geom_smooth(method = "lm",
                color = line_color,
                fill = fill_color) +
    theme_bw() +
    theme(
      legend.position = "None",
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    ) +
    xlab(xlabel) +
    ylab(ylabel)
  
  if (!is.null(corr.method)) {
    p <- p + stat_cor(method = corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}


plot_beta_div <- function(ps,
                          dist_matrix,
                          ordination,
                          group_variable,
                          add_ellipse = FALSE,
                          cols = NULL) {
  
  
  if (is.null(cols)) {
    cols <- calc_pal(ps, group_variable)
  }
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  
  # check sample data dimensions >1 otherwise phyloseq
  # fails to colour samples due to bug:https://github.com/joey711/phyloseq/issues/541
  if (dim(sample_data(ps))[2] < 2) {
    # add repeat of first column as dummy
    sample_data(ps)[, 2] <- sample_data(ps)[, 1]
  }
  plot <- plot_ordination(ps, ordination , color = group_variable)
  plot$layers[[1]] <- NULL
  
  plot_out <- plot + geom_point(size = 3, alpha = 0.75) +
    theme_bw() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    labs(caption = bquote(Adonis ~ R ^ 2 ~ .(round(ad$R2[1], 2)) ~
                            ~ p - value ~ .(ad$`Pr(>F)`[1])))
  if (add_ellipse == TRUE){
    plot_out <- plot_out + 
      geom_polygon(stat = "ellipse", aes(fill = .data [[ group_variable ]] ), alpha = 0.3)
  }
  
  return(plot_out)
}


plot_da <- function(ancom_res, groups, cols = NULL) {
  if (is.null(cols)) {
    cols <- pal_npg()(length(groups))
  }
  
  ancom_da_plot <- ancom_res %>%
    mutate(enriched_in = ifelse(Log2FC > 0, groups[2],
                                groups[1]))
  
  ancom_da_plot_sort <- ancom_da_plot %>%
    arrange(Log2FC) %>%
    mutate(Highest_classified = factor(Highest_classified, levels = Highest_classified))
  
  
  p <-
    ggplot(ancom_da_plot_sort,
           aes(x = Highest_classified, y = Log2FC, color = enriched_in)) +
    geom_point(size = 5) +
    labs(y = paste("\nLog2 Fold-Change", groups[1], "vs", groups[2], sep =
                     " "),
         x = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.position = "none"
    ) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_color_manual(values = cols)
  
  print(p)
}
