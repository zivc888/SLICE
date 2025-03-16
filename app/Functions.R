# Define helper functions
box_profile <- function(profile_mat, mut_mat, driver_gene, type = "pan", cancer_type = NA, goi = NA,
                        dist = F, DD_thresh = 0.2, highlight_cells = NA, paper_version = F){
  tmp <- profile_mat[driver_gene,] != 0
  potential_partners <- names(tmp)[tmp]
  
  # If there are no potential parnters for the driver-gene, we have nothing to show...
  if(length(potential_partners) == 0){
    return(NA)
  }
  
  # Pcell - cells which are positive for a mutation in the driver-gene (MUT).
  # Ncell- cells which are negative for a mutation in the driver-gene (WT).
  tmp <- mut_mat[driver_gene,] == 1
  Pcell <- names(tmp)[tmp]
  tmp <- mut_mat[driver_gene,] == 0
  Ncell <- names(tmp)[tmp]
  
  Pcell <- gsub("\\.", "-", Pcell)
  Ncell <- gsub("\\.", "-", Ncell)
  
  # In a cancer-specific analysis, filter cell-lines that are not relevant to us.
  if(type == "specific"){
    Pcell <- c(Pcell[METADATA[Pcell,]$primary_disease == cancer_type])
    Ncell <- c(Ncell[METADATA[Ncell,]$primary_disease == cancer_type])
  }
  
  # If there are values in the goi (gene of interest) parameter, act as if only the goi(s) is/are the potential partners to the driver-gene.
  if(length(goi) >= 1 || !is.na(goi)){
    potential_partners <- goi
  }
  
  # Iterate through all the potential partners.
  # For each - calculate the dependency of g (the current potential partner) in each cell line.
  # If there are less than 5 cell-lines with/without a mutation in the driver-gene that also have dependency data for g, return NA.
  # Else, we calculate the mean dependency of each group.
  #   *If the calculation of the means failed for some reason, we send a ggplot object containing an error message for the shiny app.
  # If there are multiple genes inside goi, or if there is 1 gene and dist is set to TRUE, we check whether g is found in goi. If it does,
  #   we create pie-charts for the composition of cancer-types in the different parts of the boxplot.
  # Regardless of the last step, we create the dependency boxplot and save/return it (goi is empty - save, goi is not empty - return).
  for(g in potential_partners){ 
    Pcell_dependencies <- CRISPR %>% filter(depmap_id %in% Pcell) %>% filter(gene_name == g)
    Pcell_dependencies <- Pcell_dependencies[,c("depmap_id", "dependency")]
    
    Ncell_dependencies <- CRISPR %>% filter(depmap_id %in% Ncell) %>% filter(gene_name == g)
    Ncell_dependencies <- Ncell_dependencies[,c("depmap_id", "dependency")]
    
    # Check for sufficient number of cell-lines in the driver-gene-MUT and the driver-gene-WT groups after filtering for a specific partner.
    if(isTRUE(nrow(Pcell_dependencies) >= 5) & isTRUE(nrow(Ncell_dependencies) >= 5)){ #5? sure?
      Pcell_dependencies$group <- "MUT"
      P_mean <- mean(Pcell_dependencies$dependency)
      
      Ncell_dependencies$group <- "WT"
      N_mean <- mean(Ncell_dependencies$dependency)
      
      # Checking if the mean calculation succeeded.
      # In the past we used to also check whether the mean-difference was higher than the DD_thresh parameter.
      if(is.na(N_mean) || is.na(P_mean) ){#|| (N_mean-P_mean < DD_thresh)){
        print("is.na(N_mean) || is.na(P_mean)")
        graph <- ggplot(data.frame(), aes(1, 1)) +
          geom_text(aes(label = "Error: not enough instances in database to create graph"),
                    size = 10, hjust = 0.5, vjust = 0.5) +
          theme_void()  # Remove axis and background elements
        
        return(graph)
      }
      
      p_num <- nrow(Pcell_dependencies)
      n_num <- nrow(Ncell_dependencies)
      data <- full_join(Pcell_dependencies, Ncell_dependencies, by=colnames(Pcell_dependencies))
      if(type == "pan"){
        data$cancer_types <- METADATA[data$depmap_id,]$primary_disease
      }
      else{
        # This is not really "cancer-types", it's "cancer-subtypes".
        data$cancer_types <- METADATA[data$depmap_id,]$lineage_subtype
      }
      
      # If there are multiple genes in the goi parameter or if there is 1 gene in it and the dist parameter is set tot TRUE.
      if((length(goi) > 1 || !is.na(goi)) && dist){
        # If the current fene is in goi, we create pie-charts for it.
        if(g %in% goi){
          pies <- list()
          pies[[1]] <- table((data %>% filter(dependency < -0.5) %>% filter(group == "WT"))$cancer_types)
          pies[[1]] <- data.frame(pies[[1]])
          
          pies[[2]] <- table((data %>% filter(dependency > -0.5) %>% filter(group == "WT"))$cancer_types)
          pies[[2]] <- data.frame(pies[[2]])
          
          pies[[3]] <- table((data %>% filter(dependency < -0.5) %>% filter(group == "MUT"))$cancer_types)
          pies[[3]] <- data.frame(pies[[3]])
          
          pies[[4]] <- table((data %>% filter(dependency > -0.5) %>% filter(group == "MUT"))$cancer_types)
          pies[[4]] <- data.frame(pies[[4]])
          
          for(ind in 1:4){
            pies[[ind]] <- ggplot(pies[[ind]], aes(x="", y=Freq, fill=Var1)) +
              geom_bar(width=1, stat="identity", color="white") +
              coord_polar("y", start=0) +
              theme_void() +
              scale_fill_manual(name = "Cancer Types",
                                values = c('#FF0000', '#A98307', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF', '#FF7D00', '#0082FF', '#7D00FF', '#FF0082', '#FF8080', '#666666', '#80FF80', '#FCFC7D','#390056', '#8080FF', '#000000', '#800000', '#808000', '#008000', '#008080', '#000080', '#800080', '#800040', '#00F6B3', '#002256', '#005655', '#564200', '#560000','#E55137','#D0D0D0'))
            plot(pies[[ind]])
          }
        }
      }
      
      data$name <- NA
      for(i in rownames(data)){
        data[i,]$name <- (METADATA %>% filter(depmap_id == data[i,]$depmap_id))$stripped_cell_line_name
      }
      
      # Creating the actual boxplot.
      DD_value <- format(round(N_mean-P_mean, 2), nsmall = 2)
      WT_percentage <- round(n_num/(p_num+n_num)*100, 2)
      MUT_percentage <-round(p_num/(p_num+n_num)*100, 2)
      box <- box_plot_creation(data, type, cancer_type, driver_gene, g, highlight_cells, paper_version,
                               DD_value, WT_percentage, MUT_percentage)
      
      # If there is no goi (we are iterating through all the potential partners), save the boxplot to a file.
      if(is.na(goi)){
        png(paste(c("/bigdata/zivcohen/graphs/boxplots_",type,"/",
                    ifelse(type=="specific",paste(c(make.names(cancer_type),": "),collapse=""),""),
                    driver_gene,"_",g,".png"), collapse=""), height = 480, width = 2*480)
        plot(box)
        dev.off()
      }
      else{
        # If there is a goi, return the boxplot.
        return(box)
      }
    }
    else{
      # This happens when there are not enough cell-lines.
      print(Pcell_dependencies)
      print(Ncell_dependencies)
      
      return(NA)
    }
  }
  
  
}

box_plot_creation <- function(data, type, cancer_type, driver_gene, partner_gene,
                              highlight_cells, paper_version, DD_value, WT_percentage, MUT_percentage){
  if(!paper_version){
    box <- ggplot(data, aes(x=group, y=dependency, fill=group)) +
      geom_boxplot(outlier.shape = NA, show.legend = TRUE) +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(aes(color=cancer_types), size=2, alpha=0.9) +
      stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
      theme_ipsum() +
      theme(plot.title = element_text(size=11),
            legend.text = element_text(size=10),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.position = "right") +
      geom_signif(comparisons = list(c("WT","MUT")),test = "t.test") +
      ggtitle(label = paste(c(partner_gene, " dependency in ", driver_gene, "+/- ",
                              ifelse(type=="specific", cancer_type, ""),
                              " cells \nDD = ", DD_value), collapse=""),
              subtitle = paste(c("WT: ", WT_percentage,
                                 "%\nMUT: ", MUT_percentage,"%"),collapse="")) +
      geom_point(data = data[data$depmap_id %in% highlight_cells,], color = "#000080") +
      geom_text(aes(label = name),
                data = data[data$depmap_id %in% highlight_cells,],
                color = "#000080", hjust = -.1)
    
    return(box)
  }
  else{
    box <- ggplot(data, aes(x=group, y=dependency, fill=group)) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      geom_jitter(aes(color=cancer_types), size=2, alpha=0.9) +
      stat_summary(fun=mean, geom="point", shape=20, size=5, color="red", fill="red") +
      theme_ipsum() +
      theme(plot.title = element_text(size=15),
            legend.text = element_text(size=10),
            axis.title.x = element_text(size = 14, hjust = 0.5),
            axis.title.y = element_text(size = 14, hjust = 0.5),
            legend.position = "none") +
      xlab(paste(c(driver_gene, "'s Mutation Status"), collapse="")) +
      ylab(paste(c(partner_gene, "'s Dependency Score"), collapse="")) +
      geom_signif(comparisons = list(c("WT","MUT")),test = "t.test") +
      ggtitle(label = paste(c(partner_gene, " dependency in ", driver_gene, "+/- ",
                              ifelse(type=="specific", cancer_type, ""),
                              " cells \nDD = ", DD_value), collapse=""),
              subtitle = "") +
      geom_point(data = data[data$depmap_id %in% highlight_cells,], color = "#000080") +
      geom_text(aes(label = name),
                data = data[data$depmap_id %in% highlight_cells,],
                color = "#000080", hjust = -.1)
    
    return(box)
  }
}
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}