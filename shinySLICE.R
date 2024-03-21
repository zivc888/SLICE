#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#load("/bigdata/zivcohen/ShinyData/shinyWorkspace.RData")

library(shiny)
library(dplyr)
library(shinycssloaders)
library(plotly)
library(network)
library(GGally)
library(sna)
library(dichromat)
library(bslib)
library(viridis)
library(hrbrthemes)
library(ggpubr)
library(reshape2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bs_theme(bootswatch = bootswatch_themes()[7]),
  
  # Application title
  titlePanel(
    "SLICE - Synthetic Lethality Identification through CRISPR-based Evaluation"
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "analysisType",
        "Analysis Type",
        c(
          "Choose",
          "Mutation-based Analysis",
          "Expression-based Analysis"
        )
      ),
      selectInput(
        "analysisScope",
        "Analysis Scope",
        c("Choose", "Pan-cancer Analysis", "Cancer-specific Analysis")
      ),
      conditionalPanel(
        condition = "input.analysisScope == 'Cancer-specific Analysis'",
        selectInput(
          "cancerType",
          "Select Cancer Type",
          choices = sort(cancer_types_filt)
        )
      ),
      
      conditionalPanel(
        condition = "input.analysisType == 'Mutation-based Analysis'",
        selectInput("driverGene", "Select Driver Gene", choices = sort(common_genes))
      ),
      
      conditionalPanel(
        condition = "input.analysisType == 'Expression-based Analysis'",
        sliderInput(
          "corrThreshold",
          "Correlation Threshold",
          0.5,
          max = 1,
          min = 0,
          step = 0.05
        )
      ),
      
      textInput(
        "highlight",
        "Highlight by gene name",
        value = "",
        width = NULL,
        placeholder = NULL
      ),
      checkboxInput("match", "exact match?"),
      checkboxInput("onlydrug", "Show Only Druggable Genes."),
      
      actionButton("plot1", "Plot It!"),
      
      conditionalPanel(
        condition = "input.analysisType == 'Expression-based Analysis'",
        checkboxInput("showpath", "Show the path between the following nodes:"),
      ),
      
      conditionalPanel(
        condition = "input.showpath",
        selectInput(
          "path1",
          "",
          choices = "None"
        ),
        selectInput(
          "path2",
          "",
          choices = "None"
        )
      )
    ),
    
    mainPanel(
      plotlyOutput("generalPlot") %>% withSpinner(color = "#0066aa"),
      plotOutput("specificPlot") %>% withSpinner(color = "#0066aa")
    )
  )
)

lfc <- NULL
pvalue <- NULL

profile_matrix <- dependencyProfile_matrix
mut_matrix <- mutation_matrix

# Define server logic
server <- function(input, output, session) {
  values <- reactiveValues(
    GeneOfInterest = NULL,
    SetOfInterest = NULL,
    OnlyDruggable = FALSE,
    selectedDriverGene = NULL,
    selectedCancerType = NULL,
    analysisScope = NULL,
    analysisType = "M"
  )
  
  observe({
    # Update choices for "cancerType" based on the selected analysisScope
    if (input$analysisScope == "Cancer-specific Analysis") {
      updateSelectInput(session,
                        "driverGene",
                        choices = sort(row.names(mut_spec_lfcs[[input$cancerType]])))
    }
    else{
      updateSelectInput(session,
                        "driverGene",
                        choices = sort(common_genes))
    }
  })
  observe({
    # Update choices for "cancerType" based on the selected analysisScope
    if (input$analysisType == "Expression-based Analysis") {
      updateSelectInput(session,
                        "cancerType",
                        choices = sort(gsub("/", " ", filtered_primary_disease)))
      updateTextInput(session,
                      "highlight",
                      label = "Filter by node name")
      updateCheckboxInput(session,
                          "onlydrug",
                          label = "Show only druggable genes")
    }
    else{
      updateSelectInput(session,
                        "cancerType",
                        choices = sort(gsub("/", " ", cancer_types_filt)))
      updateTextInput(session,
                      "highlight",
                      label = "Highlight by gene name")
      updateCheckboxInput(session,
                          "onlydrug",
                          label = "Highlight druggable genes")
    }
  })
  observe({
    values$OnlyDruggable <- input$onlydrug
  })
  
  observeEvent(input$plot1, {
    values$GeneOfInterest <- NULL
    
    # Perform analysis based on analysis type
    if (input$analysisType == "Mutation-based Analysis") {
      values$analysisType <- "M"
      values$selectedDriverGene <- input$driverGene
      
      if (input$analysisScope == "Pan-cancer Analysis") {
        values$selectedCancerType <- NULL
        lfc <- lfc_matrix[values$selectedDriverGene, ]
        pvalue <- pvalue_matrix[values$selectedDriverGene, ]
      }
      else if (input$analysisScope == "Cancer-specific Analysis") {
        values$selectedCancerType <- input$cancerType
        lfc <-
          mut_spec_lfcs[[values$selectedCancerType]][values$selectedDriverGene, ]
        pvalue <-
          mut_spec_pvalues[[values$selectedCancerType]][values$selectedDriverGene, ]
      }
    }
    else if (input$analysisType == "Expression-based Analysis") {
      values$analysisType <- "E"
      
      if (input$analysisScope == "Pan-cancer Analysis") {
        values$selectedCancerType <- NULL
      }
      else if (input$analysisScope == "Cancer-specific Analysis") {
        values$selectedCancerType <- input$cancerType
      }
    }
    
    output$generalPlot <- renderPlotly({
      if (values$analysisType == 'M') {
        df <-
          data.frame(
            "gene_name" = names(lfc),
            "LFC" = lfc,
            "qvalue" = pvalue
          )
        
        # Create the plotly object directly using plot_ly()
        vline <- function(x = 0, color = "red") {
          list(
            type = "line",
            y0 = 0,
            y1 = 1,
            yref = "paper",
            x0 = x,
            x1 = x,
            line = list(color = color, dash = "dot")
          )
        }
        hline <- function(y = 0, color = "red") {
          list(
            type = "line",
            y0 = y,
            y1 = y,
            xref = "paper",
            x0 = 0,
            x1 = 1,
            line = list(color = color, dash = "dot")
          )
        }
        
        point_color <- function(lfc, qvalue, lfc_threshold, qvalue_threshold, gene) {
          if (values$OnlyDruggable && !(gene %in% druggable_genes)) {
            return("lightgray")
          }
          
          if (abs(lfc) >= lfc_threshold && qvalue <= qvalue_threshold) {
            return("red")
          } else if (abs(lfc) >= lfc_threshold || qvalue <= qvalue_threshold) {
            return("blue")
          } else {
            return("gray")
          }
        }
        
        df$point_color <- mapply(point_color, df$LFC, df$qvalue, 0.2, 0.05, df$gene_name)
        
        volcano_plotly <- plot_ly(source = "generalPlot") %>%
          add_trace(
            data = df,
            x = ~ -LFC,
            y = ~ -log10(qvalue),
            text = ~ gene_name,
            type = "scatter",
            mode = "markers",
            key = ~ gene_name,
            marker = list(
              size = 5,
              color = ~ point_color,
              opacity = 0.6
            )
          ) %>%
          layout(
            shapes = list(
              vline(0.2, color = "green"),
              vline(-0.2, color = "green"),
              hline(-log10(0.05), color = "green")
            ),
            title = paste(
              values$selectedDriverGene,
              "-WT vs. ",
              values$selectedDriverGene,
              "-mutated",
              if (is.null(values$selectedCancerType))
                ""
              else
                paste(" in ", values$selectedCancerType, " cells", sep = ""),
              sep = ""
            ),
            xaxis = list(title = "Delta-Dependency (DD)"),
            yaxis = list(title = "-log10(p-value)"),
            hoverlabel = list(
              bgcolor = "white",
              font = list(family = "Arial, sans-serif")
            )
          )
        
        if (input$highlight != "") {
          if(!input$match){
            highlight_df <-
              df[grepl(paste(strsplit(input$highlight, split = ", ")[[1]],
                             collapse = '|'),
                       df$gene_name,
                       ignore.case = TRUE), ]
          }
          else{
            highlight_df <-
              df[df$gene_name %in% strsplit(input$highlight, split = ", ")[[1]],]
          }
          
          if (length(row.names(highlight_df)) != 0) {
            volcano_plotly <-
              volcano_plotly %>% add_annotations(
                x = -highlight_df$LFC,
                y = -log10(highlight_df$qvalue),
                text = row.names(highlight_df),
                xref = "x",
                yref = "y",
                showarrow = TRUE,
                arrowhead = 4,
                arrowsize = .5,
                standoff = 10,
                textangle = -60,
                ax = 20,
                ay = -40
              )
          }
        }
        
        return(volcano_plotly %>% event_register("plotly_click"))
      }
      else{
        # Define the color mapping function
        custom_color_map <- function(y) {
          a <- 1 / (1 + exp(20 * (abs(y) - 0.6)))
          
          if (y > 0) {
            return(rgb(1, 0, 0, alpha = 1 - a))
          }
          if (y < 0) {
            return(rgb(0, 0, 1, alpha = 1 - a))
          }
          else{
            return("white")
          }
        }
        
        # Create a data frame from the filtered matrix
        edge_df <- correlation_dfs[[ifelse(!is.null(values$selectedCancerType),
                                           values$selectedCancerType,
                                           "pan")]]
        edge_df <- edge_df[, 1:3]
        colnames(edge_df) <- c("gene_set", "gene_name", "corr")
        edge_df[abs(edge_df$corr) < input$corrThreshold, ]$corr <- 0
        
        all_edges <- NULL
        if (input$highlight != "") {
          highlights <- strsplit(input$highlight, split = ", ")[[1]]
          grep_highlights <- strsplit(input$highlight, split = ", ")[[1]]
          
          if(!input$match){
            highlight_df <- rbind(edge_df[grepl(paste(highlights,collapse="|"),
                                                edge_df$gene_name,
                                                ignore.case = TRUE), ],
                                  edge_df[grepl(paste(highlights,collapse="|"),
                                                edge_df$gene_set,
                                                ignore.case = TRUE), ])
          }
          else{
            highlight_df <- rbind(edge_df[edge_df$gene_name %in% highlights, ],
                                  edge_df[edge_df$gene_set %in% highlights, ])
          }
          
          if (length(row.names(highlight_df)) != 0) {
            all_edges <- edge_df
            edge_df <- highlight_df
          }
        }
        else{
          all_edges <- edge_df
        }
        
        if(input$showpath){
          relevant_nodes <- unique(c(sort(edge_df[edge_df$corr != 0,]$gene_set),
                                     sort(edge_df[edge_df$corr != 0,]$gene_name)))
          
          if(input$path1 == "None" && input$path2 == "None"){
            updateSelectInput(session,
                              "path1",
                              choices = c("None", relevant_nodes)
            )
            updateSelectInput(session,
                              "path2",
                              choices = c("None", relevant_nodes)
            )
          }
          
          
          if(input$path1 != "None" && input$path2 != "None"){
            highlight_nodes <- c(input$path1, input$path2)
            print(highlight_nodes)
            all_possible_nodes <- subset(all_edges, corr != 0)
            
            # TAKE INTO CONSIDERATION THE WEIGHTS?
            FindShortestPath <-
              function(all_possible_nodes, source, end, path, depth, thresh){
                if(source %in% all_possible_nodes$gene_name){
                  Neighbors <-
                    unique((all_possible_nodes %>% filter(gene_name == source))$gene_set)
                }
                else if(source %in% all_possible_nodes$gene_set){
                  Neighbors <-
                    unique((all_possible_nodes %>% filter(gene_set == source))$gene_name)
                }
                else{
                  print("hmmm...")
                  #oops
                }
                
                Neighbors <- Neighbors[!(Neighbors %in% path)]
                
                if(depth == thresh){
                  if(!(end %in% Neighbors)){
                    return(NULL)
                  }
                  else{
                    return(c(path, list(end)))
                    
                  }
                }
                
                for(n in Neighbors){
                  newpath <-
                    FindShortestPath(all_possible_nodes,
                                     n,
                                     end,
                                     c(path, list(n)),
                                     depth+1,
                                     thresh)
                  if(!is.null(newpath)){
                    return(newpath)
                  }
                }
                
                return(NULL)
              }
            
            
            if(length(highlight_nodes) == 2){
              
              path <- NULL
              t <- 0
              while(is.null(path) && t <= 5){
                path <- FindShortestPath(all_possible_nodes,
                                         highlight_nodes[1],
                                         highlight_nodes[2],
                                         list(highlight_nodes[1]),
                                         0,
                                         t)
                t <- t+1
              }
              
              if(!is.null(path)){
                nodes <- unlist(path)
                print(nodes)
                
                for(i in 2:(length(path)-1)){
                  if(nodes[i] %in% edge_df$gene_name){
                    edge_df <- rbind(edge_df,
                                     all_possible_nodes %>%
                                       filter(gene_name == nodes[i]) %>%
                                       filter(gene_set == nodes[i+1]))
                  }
                  else{
                    edge_df <- rbind(edge_df,
                                     all_possible_nodes %>%
                                       filter(gene_set == nodes[i]) %>%
                                       filter(gene_name == nodes[i+1]))
                  }
                  
                }
              }
            }
          }
        }
        
        if (input$onlydrug) {
          edge_df <- edge_df %>% filter(gene_name %in%
                                          intersect(edge_df$gene_name, druggable_genes))
        }
        
        # Filter out rows with zero values (no edges)
        filtered_edges <- subset(edge_df, corr != 0)
        filtered_edges <- unique(filtered_edges)
        
        if (length(row.names(filtered_edges)) == 0) {
          graph <- ggplot(data.frame(), aes(1, 1)) +
            geom_text(
              aes(label = "Error: no cell fitted the filter, try using a lower threshold"),
              size = 5,
              hjust = 0.5,
              vjust = 0.5
            ) +
            theme_void()  # Remove axis and background elements
          
          return(ggplotly(graph) %>% event_register("plotly_click"))
        }
        
        # Calculate the color values based on the correlation values
        filtered_edges$color <-
          sapply(filtered_edges$corr, custom_color_map)
        
        # Get unique node names
        all_node_names <-
          unique(c(filtered_edges$gene_name, filtered_edges$gene_set))
        
        all_node_names <-
          unique(c(filtered_edges$gene_name, filtered_edges$gene_set))
        original_matrix <-
          funrar::stack_to_matrix(filtered_edges, "gene_set", "gene_name", "corr")
        
        adjacency_matrix <-
          matrix(
            0,
            nrow = length(all_node_names),
            ncol = length(all_node_names),
            dimnames = list(all_node_names, all_node_names)
          )
        
        for (i in 1:nrow(original_matrix)) {
          for (j in 1:ncol(original_matrix)) {
            if (!is.na(original_matrix[i, j])) {
              adjacency_matrix[row.names(original_matrix)[i], colnames(original_matrix)[j]] <-
                original_matrix[i, j]
            }
          }
        }
        
        netmat1_matrix <- adjacency_matrix
        
        # Remove NA values from the matrix
        df <- na.omit(netmat1_matrix)
        
        # Create a network object from the adjacency matrix
        net1 <-
          network(df, matrix.type = "adjacency", loops = FALSE)
        
        # Create a network plot using ggnet2
        network <-
          ggnet2(
            net1,
            node.size = 5,
            node.color = rep(c("orange", "purple"),
                             times = c(
                               length(unique(filtered_edges$gene_name)),
                               length(unique(filtered_edges$gene_set))
                             )),
            edge.size = 1,
            edge.color = filtered_edges$color,
            label = all_node_names
          )
        
        # Convert the ggplot plot to a Plotly plot
        plotly_network <-
          ggplotly(network, source = "generalPlot", tooltip = c())
        
        plotly_network <- plotly_network %>% add_trace(
          x = network$data$x,
          y = network$data$y,
          hoverinfo = "text",
          text = all_node_names,
          key = all_node_names,
          type = "scatter",
          mode = "markers",
          marker = list(
            size = 25,
            color = rep(c("orange", "purple"),
                        times = c(
                          length(unique(filtered_edges$gene_name)),
                          length(unique(filtered_edges$gene_set))
                        )),
            opacity = 0
          )
        )
        
        # Customize layout
        layout <- plotly_network %>%
          layout(
            title = paste(
              "Combinations of genes and gene-sets",
              if (is.null(values$selectedCancerType))
                ""
              else
                paste(" in ", values$selectedCancerType, " cells", sep = ""),
              sep = ""
            ),
            font = list(size = 12),
            paper_bgcolor = "white"
          )
        
        return(layout %>% event_register("plotly_click"))
      }
    })
    
  })
  
  # Respond to clicking the main plot
  observeEvent(event_data(event = "plotly_click",
                          source = "generalPlot"), {
                            clicked <- event_data(event = "plotly_click",
                                                  source = "generalPlot")
                            if (!is.null(clicked)) {
                              if (values$analysisType == "E" &&
                                  (clicked$key %in% names(gene_sets))) {
                                values$SetOfInterest <- clicked$key
                              }
                              else{
                                values$GeneOfInterest <- clicked$key
                              }
                            }
                          })
  
  output$specificPlot <- renderPlot({
    paper_version <- T
    
    if (values$analysisType == "M") {
      validate(need(values$GeneOfInterest, message = "click on a dot"))
      
      mut_matrix <- mutation_matrix
      
      if (is.null(values$selectedCancerType)) {
        # pan
        profile_matrix <- dependencyProfile_matrix
        values$analysisScope <- "pan"
        currentCancerType <- NULL
      }
      else{
        # specific
        profile_matrix <-
          mut_spec_dependencyProfiles[[values$selectedCancerType]]
        values$analysisScope <- "specific"
        
        if (values$selectedCancerType == "Colon Colorectal Cancer") {
          currentCancerType <- "Colon/Colorectal Cancer"
        }
        else{
          currentCancerType <- values$selectedCancerType
        }
      }
      
      print(values$GeneOfInterest)
      
      box <- box_profile(
        profile_mat = profile_matrix,
        mut_mat = mut_matrix,
        type = values$analysisScope,
        cancer_type = currentCancerType,
        driver_gene = values$selectedDriverGene,
        goi = values$GeneOfInterest,
        DD_thresh = 0,
        paper_version = paper_version
      )
      # highlight_cells = (METADATA %>% filter(cell_line_name %in% c("5637", "T24")))$depmap_id
      
      graph <- box
    }
    else{
      validate(need(values$GeneOfInterest, message = "click on a gene"))
      validate(need(values$SetOfInterest, message = "click on a set"))
      
      print(values$GeneOfInterest)
      print(values$SetOfInterest)
      
      gene <- values$GeneOfInterest
      gene_set <- values$SetOfInterest
      
      if(is.null(values$selectedCancerType)){
        current_depmap_ids <- METADATA$depmap_id
      }
      else{
        if (values$selectedCancerType == "Colon Colorectal Cancer") {
          currentCancerType <- "Colon/Colorectal Cancer"
        }
        else{
          currentCancerType <- values$selectedCancerType
        }
        
        current_depmap_ids <-
          METADATA$depmap_id[METADATA$primary_disease == currentCancerType]
      }
      
      tmp <-
        TPM[, c("depmap_id", "cell_line", "gene_name", "rna_expression")] %>%
        filter(depmap_id %in% current_depmap_ids)
      
      tmp <-
        tmp[tmp[, "gene_name"] == gene, ][, c("depmap_id", "cell_line", "rna_expression")]
      
      tmp1 <-
        CRISPR[CRISPR$gene_name == gene, ][, c("depmap_id", "dependency", "gene_name")]
      tmp <- merge(tmp,
                   tmp1[tmp1$gene_name == gene, ][, c("depmap_id", "dependency")],
                   by = "depmap_id")
      
      tmp1 <-
        MUT[MUT$gene_name == gene, ][, c("depmap_id", "var_annotation", "gene_name")]
      tmp <- left_join(tmp,
                       tmp1[tmp1$gene_name == gene, ][, c("depmap_id", "var_annotation")],
                       by = "depmap_id")
      
      tmp1 <-
        METADATA[METADATA$depmap_id %in% tmp$depmap_id, ][, c("depmap_id",
                                                              "lineage_subtype")]
      tmp <- left_join(tmp,
                       tmp1,
                       by = "depmap_id")
      
      tmp[is.na(tmp$var_annotation), ]$var_annotation <- "unknown"
      
      tmp1 <-
        funrar::stack_to_matrix(
          TPM %>% filter(depmap_id %in% current_depmap_ids),
          'gene_name',
          'cell_line',
          'rna_expression'
        )
      ssgsea_scores <- ssgsea(tmp1, gene_sets)
      enrichment_df <- melt(
        ssgsea_scores,
        varnames = c("gene_set", "cell_line"),
        value.name = "Enrichment"
      )
      
      tmp <-
        left_join(tmp, enrichment_df[enrichment_df$gene_set == gene_set, ],
                  by = c("cell_line" = "cell_line"))
      
      if(paper_version){
        graph <-
          ggplot(tmp,
                 aes(x = Enrichment, y = dependency))
      }
      else{
        graph <-
          ggplot(tmp,
                 aes(x = Enrichment, y = dependency, col = var_annotation))
      }
      graph <- graph +
        theme_minimal() +
        geom_point() +
        geom_smooth(
          method = glm ,
          color = "red",
          fill = "#69b3a2",
          se = TRUE
        ) +
        scale_fill_viridis(option = "inferno") +
        ggtitle(
          label = ifelse(paper_version, "", ifelse(!is.null(values$selectedCancerType),
                                                   values$selectedCancerType,
                                                   "Pan")),
          subtitle = paste(gene, "'s dependency and ", gene_set, "'s enrichment", sep = "")
        ) +
        stat_cor(method = "pearson")
      
    }
    
    ggsave(plot = graph, filename = paste("/bigdata/zivcohen/graphs/promising_pSLPs/",
                                          ifelse(values$analysisType == "M", "BOX: ", "REGRESSION: "),
                                          values$selectedCancerType,"_",
                                          ifelse(values$analysisType == "M",
                                                 paste(values$selectedDriverGene, "_", values$GeneOfInterest, sep=""),
                                                 paste(values$GeneOfInterest, "_", values$SetOfInterest, sep="")),
                                          ".png", sep=""),
           height = ifelse(values$analysisType == "M",3000,1500), width = 1500, limitsize = FALSE, units = "px")
    return(graph)
  })
}


# Run the application
shinyApp(ui = ui, server = server)
