# Required Packages 
library(shiny)
library(dplyr)
library(tidyr)
library(shinyWidgets)
library(shinythemes)
library(markdown)
library(vroom)
library(ggpubr)
library(ggplot2)
library(readr)
library(DT)


# Loading necessary data for app functionality # 
# Load in the required tables for HUMAN HOST NESTED GENES
hn_genes <- vroom("hn_human_master.csv")
# Addition of column detailing coordinates
hn_genes$x <- paste0(hn_genes$geneChr_Host, ":", hn_genes$geneStart_Nested)
hn_genes$Position_Nested <- paste0(hn_genes$x , "-", hn_genes$geneEnd_Nested)
hn_genes$x <- paste0(hn_genes$geneChr_Host, ":", hn_genes$geneStart_Host)
hn_genes$Position_Host <- paste0(hn_genes$x , "-", hn_genes$geneEnd_Host)
hn_genes <- hn_genes %>% select(-x)
# Selection / Filtering of table to show desired information on the app 
hn_genes_selected <- hn_genes %>% 
  select(Position_Nested, ENST_Nested, SYMBOL_Nested, biotype_Nested, 
         Position_Host, ENST_Host, SYMBOL_Host, biotype_Host,
         pairid) %>% 
  dplyr::rename(pair_id = pairid)
# Loading and sorting expression data 
expression_data <- vroom("expression_data.csv", col_select = c(-1))%>% mutate(nest = log(nest+1)) %>% mutate(host=log(host+1))
expression_data2 <- 
  expression_data %>% gather("nest_host", "lognorm", -sample, -pair_id, -tissue ) %>% 
  left_join(hn_genes_selected) %>% mutate(GeneID = case_when(
    endsWith(nest_host, "nest") ~ .$SYMBOL_Nested,
    endsWith(nest_host, "host") ~ .$SYMBOL_Host
  )) %>% select(2:5,14)


# Load in the required tables for MOUSE DATA, performed in similar manner to above
mouse_hn_genes <- vroom("hn_mouse_master.csv")
mouse_hn_genes$x <- paste0(mouse_hn_genes$geneChr_Host, ":", mouse_hn_genes$geneStart_Nested)
mouse_hn_genes$Position_Nested <- paste0(mouse_hn_genes$x , "-", mouse_hn_genes$geneEnd_Nested)
mouse_hn_genes$x <- paste0(mouse_hn_genes$geneChr_Host, ":", mouse_hn_genes$geneStart_Host)
mouse_hn_genes$Position_Host <- paste0(mouse_hn_genes$x , "-", mouse_hn_genes$geneEnd_Host)
mouse_hn_genes <- mouse_hn_genes %>% select(-x)
mouse_hn_genes_selected <- mouse_hn_genes %>% 
  select(Position_Nested, ENST_Nested, SYMBOL_Nested, biotype_Nested, 
         Position_Host, ENST_Host, SYMBOL_Host, biotype_Host,
         pairid) %>% 
  dplyr::rename(pair_id = pairid)
mouse_expression_data <- vroom("mouse_expression_data.csv", col_select = c(-1)) %>% mutate(nest = log(nest+1)) %>% mutate(host=log(host+1))
m_expression_data2 <- mouse_expression_data %>% gather("nest_host", "lognorm", -sample, -pair_id, -tissue ) %>% 
  left_join(mouse_hn_genes_selected) %>% mutate(GeneID = case_when(
    endsWith(nest_host, "nest") ~ .$SYMBOL_Nested,
    endsWith(nest_host, "host") ~ .$SYMBOL_Host
  )) %>% select(2:5,14)

# RNA-Sequencing Metadata Table
hg19_meta <- read_csv("ENCODE_human_metadata.csv")
mm10_meta <- read_csv("ENCODE_mouse_metadata.csv")


# Load in testis data for scRNA-seq
proportion <- vroom("All_proportions_new.csv")
# Order
proportion$Cluster_identity <- factor(proportion$Cluster_identity , levels = c("Endothelial cells","Fibrotic peritubular myoid cells","Leydig cells",
                                                                               "Macrophages","PMCs","Sertoli cells",
                                                                               "Undifferentiated spermatogonia","Diff.spermatogonia/Preleptotene",
                                                                               "Leptotene","Zygotene","Pachytene",
                                                                               "Diplotene","Meiotic divisions",
                                                                               "Early spermatids","Late spermatids"))
# Filtering 
proportion2 <-  proportion %>% select(-cluster, -not_expressing_prop) %>%
  pivot_longer(cols=c("host_only_prop","nested_only_prop","co_exp_prop"),names_to = "type", values_to = "value")
# Obtain pair_id information
hg19_hngenes <- vroom("hn_human_master.csv") %>% select(SYMBOL_Host, SYMBOL_Nested, pairid) %>% 
  mutate(H_N = paste(SYMBOL_Host, ":", SYMBOL_Nested))
# Generate proportion table 
proportion3 <- merge(proportion2, hg19_hngenes, by = "pairid")

# Loading in testis expression data 
Testis_nested_data <- as.data.frame(vroom("testis_nested_data.csv"))
rownames(Testis_nested_data) <- as.vector(Testis_nested_data[,1])
Testis_nested_data <- Testis_nested_data[,-1]

Testis_host_data <- as.data.frame(vroom("testis_host_data.csv"))
rownames(Testis_host_data) <- as.vector(Testis_host_data[,1])
Testis_host_data <- Testis_host_data[,-1]

tissue_meta <- vroom("tissue_meta.csv")

pair_meta <- proportion3 %>% select(pair_id=pairid , SYMBOL_Host, SYMBOL_Nested, H_N)



shinyServer(function(input, output, session) {
  
  ### ### ### ### ### ### ### ### 
  # Human Data #
  ### ### ### ### ### ### ### ### 
  
  # Viewing and selecting Host Nested Genes from the table # 
  table <- reactive(hn_genes_selected)
  selection <- reactive(table()[input$table_rows_selected,] %>% pull(pair_id))
  host_gene <- reactive(table()[input$table_rows_selected,] %>% pull(SYMBOL_Host))
  nested_gene <- reactive(table()[input$table_rows_selected,] %>% pull(SYMBOL_Nested))
  
  output$table  <-  DT::renderDataTable( table(), selection = 'single')
  
  output$selection_text <- renderText({
    paste("The Host / Nested Gene pair you have selected is: ", "\n", 
          "Host Gene: ", host_gene(), "\n", 
          "Nested Gene: ", nested_gene())
  })
  
  
  
  # Expression of genes across tissues # 
  tissue_plot <- reactive(
    expression_data2 %>% filter(pair_id %in% selection()) %>% 
      ggplot(aes(tissue, lognorm, fill = tissue)) + 
      geom_bar(stat = "summary") + ylab("log_normalised_counts") + xlab("") + 
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      facet_grid(GeneID ~.)
  )
  
  correlation_plot <- reactive(expression_data %>% filter(pair_id %in% selection()) %>% 
                                 ggscatter(x = "nest", y = "host", fill = "tissue", color = "tissue", 
                                           add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
                                 stat_cor(method = "spearman", label.y.npc="top", label.x.npc = "left") + 
                                 theme(legend.position="none") + xlab(nested_gene()) + ylab(host_gene())
  )
  
  
  plot1 <- eventReactive(input$plot_button,{
    correlation_plot() 
  }
  )
  
  plot2 <- eventReactive(input$plot_button2,{
    tissue_plot()
  }
  )
  
  output$expression <- renderPlot({
    plot1()
  }, res = 96, height = 300, width = 300)
  
  output$expression2 <- renderPlot({
    plot2()
  }, res = 96, height = 300, width = 500)
  
  
  
  
  # Generation of Custom UCSC link based on genes # 
  
  # Attempting with fixed link, where you specify your session in the URL, hgsid changes depending on user
  hg19_chr <- reactive(hn_genes %>% filter(pairid == selection()) %>% pull(geneChr_Host))
  hg19_start <- reactive(hn_genes %>% filter(pairid == selection()) %>% pull(geneStart_Host))
  hg19_end <- reactive(hn_genes %>% filter(pairid == selection()) %>% pull(geneEnd_Host))
  hg19_url <- reactive(
    paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jcain-oakey&hgS_otherUserSessionName=hg19_withcomprehensive_polyA&position=",
           hg19_chr(),
           "%3A",
           hg19_start(),
           "%2D",
           hg19_end()
    )
  )
    
    url <- reactive(a("hg19 UCSC Genome Browser Link","for Host/Nested pair:",host_gene()," / ",nested_gene(), href=hg19_url()))
    url_click <- eventReactive(input$url_button, {
      url()
    })
    output$tab <- renderUI({
      tagList(url_click())
    })
    
    
    
  
  
  
  ### ### ### ### ### ### ### ### 
  # Mouse Data #
  ### ### ### ### ### ### ### ### 
  
  # Table # 
  m_table <- reactive(mouse_hn_genes_selected)
  m_selection <- reactive(m_table()[input$m_table_rows_selected,] %>% pull(pair_id))
  m_host_gene <- reactive(m_table()[input$m_table_rows_selected,] %>% pull(SYMBOL_Host))
  m_nested_gene <- reactive(m_table()[input$m_table_rows_selected,] %>% pull(SYMBOL_Nested))
  
  output$m_table  <-  DT::renderDataTable( m_table(), selection = 'single')
  
  
  output$m_selection_text <- renderText({
    paste("The Host / Nested Gene pair you have selected is: ", "\n", 
          "Host Gene: ", m_host_gene(), "\n", 
          "Nested Gene: ", m_nested_gene())
  })
  
  
  # Expression # 
  m_tissue_plot <- reactive(
    m_expression_data2 %>% filter(pair_id %in% m_selection()) %>% 
      ggplot(aes(tissue, lognorm, fill = tissue)) + 
      geom_bar(stat = "summary") + ylab("log_normalised_counts") + xlab("") + 
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      facet_grid(GeneID ~.)
  )
  
  
  
  m_correlation_plot <- reactive(mouse_expression_data %>% filter(pair_id %in% m_selection()) %>% 
                                   ggscatter(x = "nest", y = "host", fill = "tissue", color = "tissue", 
                                             add = "reg.line", add.params = list(color = "blue", fill = "lightgray")) + 
                                   stat_cor(method = "spearman", label.y.npc="top", label.x.npc = "left") + 
                                   theme(legend.position="none") + xlab(m_nested_gene()) + ylab(m_host_gene())
  )
  
  
  
  
  m_plot1 <- eventReactive(input$m_plot_button,{
    m_correlation_plot() 
  }
  )
  
  m_plot2 <- eventReactive(input$m_plot_button2,{
    m_tissue_plot()
  }
  )
  
  output$m_expression <- renderPlot({
    m_plot1()
  }, res = 96, height = 300, width = 300)
  
  output$m_expression2 <- renderPlot({
    m_plot2()
  }, res = 96, height = 300, width = 500)
  
  
  # UCSC URL # 
  mm10_chr <- reactive(mouse_hn_genes %>% filter(pairid == m_selection()) %>% pull(geneChr_Host))
  mm10_start <- reactive(mouse_hn_genes %>% filter(pairid == m_selection()) %>% pull(geneStart_Host))
  mm10_end <- reactive(mouse_hn_genes %>% filter(pairid == m_selection()) %>% pull(geneEnd_Host))
  mm10_url <- reactive(
    paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jcain-oakey&hgS_otherUserSessionName=mm10_withcomprehensive_polyA&position=",
           mm10_chr(),
           "%3A",
           mm10_start(),
           "%2D",
           mm10_end()
    )
  )
  
  
  m_url <- reactive(a("mm10 UCSC Genome Browser Link","for Host/Nested pair:",m_host_gene()," / ",m_nested_gene(), href=mm10_url()))
  m_url_click <- eventReactive(input$m_url_button, {
    m_url()
  })
  output$m_tab <- renderUI({
    tagList(m_url_click())
  })
  
    

    ### ### ### ### ### ### ### ### 
    # Spermatogenesis Tab #
    ### ### ### ### ### ### ### ### 
    
        
    df <- reactive({ 
      df <- proportion3 %>% filter(H_N == input$variable)
    })
    
    output$plot <- renderPlot({ 
      dfb <- df()
      ggplot(dfb, aes(x=Cluster_identity, y=value, group=type)) +
        geom_line(aes(color=type))+
        geom_point(aes(color=type))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_color_manual(values = c("#d9227eff","#f2a800ff","#002ca8ff")) + 
        xlab("") + ylab("Percentage of cells")
      
    })
    
    
    # Scatter Plot Facet 
    testis_pair <- reactive({ 
      pairid_temp <- proportion3 %>% filter(H_N == input$variable) %>% select(pairid) %>% pull(pairid)
      testis_pair <- pairid_temp[1]
      
    })
    
    host_id <- reactive({ 
      host_id_temp <- proportion3 %>% filter(H_N == input$variable) %>% select(SYMBOL_Host) %>% pull(SYMBOL_Host)
      host_id <- host_id_temp[1]
    })
    
    nested_id <- reactive({ 
      nested_id_temp <- proportion3 %>% filter(H_N == input$variable) %>% select(SYMBOL_Nested) %>% pull(SYMBOL_Nested)
      nested_id <- nested_id_temp[1]
      
    })
    
    
    output$plot2 <- renderPlot({ 
      testis_pair <- testis_pair()
      nested_id <- nested_id()
      host_id <- host_id()
      group1 <- Testis_nested_data[testis_pair,]
      group2 <- Testis_host_data[testis_pair,]
      data <- data.frame(nest=as.numeric(group1), host=as.numeric(group2), sample = colnames(Testis_nested_data), pair_id = testis_pair)
      data <- merge(data, tissue_meta, by = "sample")
      data <- merge(data, unique(pair_meta), by = "pair_id")
      halfmax_host <- max(as.numeric(data$host), na.rm = T)/4
      halfmax_nested <- max(as.numeric(data$nest), na.rm = T)/4
      ggplot(data, aes(nest,host)) +
        facet_wrap(.~tissue) +
        annotate("rect", xmin = halfmax_nested, xmax = Inf, ymin = halfmax_host , ymax = Inf,fill = "#f9d2e4")+
        annotate("rect",xmin = -Inf, xmax = halfmax_nested, ymin = halfmax_host , ymax = Inf,fill = "#fff0cc")+
        annotate("rect",xmin = halfmax_nested, xmax = Inf, ymin = -Inf , ymax = halfmax_host,fill = "#ccdaff")+
        geom_point(aes(colour = factor(tissue)),size = 2, alpha = 0.5) +
        geom_hline(yintercept=halfmax_host, linetype="dashed", color = "grey",  linewidth =1.5) +
        geom_vline(xintercept = halfmax_nested, linetype="dashed", color = "grey", linewidth=1.5) +
        theme(legend.position="none")+
        xlab(paste(nested_id, "expression \n log(normalised counts)")) + ylab(paste(host_id, "expression \n log(normalised counts)")) + theme_bw() + labs(color='Tissue') 
      
      
    }  )
    
    
    
  
  
  # Outputs for the About Page # 
  
  output$pipeline <- renderImage({
    filename <- "hn_gene_pipeline.png"
    list(src = filename, width = "1000", height = "450")
  }, deleteFile = FALSE)
  
  
  output$downloadMouseMetaData <- downloadHandler(
    filename = function() { 
      paste("mouse_ENCODE_metadata_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(mm10_meta, file)
    })
  output$downloadHumanMetaData <- downloadHandler(
    filename = function() { 
      paste("human_ENCODE_metadata_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(hg19_meta, file)
    })
  
  
  output$downloadMouseHnGenes <- downloadHandler(
    filename = function() { 
      paste("mm10_hostnested_genes", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(mouse_hn_genes, file)
    })
  output$downloadHumanHnGenes <- downloadHandler(
    filename = function() { 
      paste("hg19_hostnested_genes", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(hn_genes, file)
    })
  

  output$UCSC <- renderImage({
    filename <- "UCSC_custom_track_explanation.png"
    list(src = filename, width = "1300", height = "500")
  }, deleteFile = FALSE)
  
}
)


paste0("nested_id","expression \n log(normalised counts)")
paste("nested_id", "expression \n log(normalised counts)")
