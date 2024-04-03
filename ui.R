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




shinyUI(navbarPage(title = div(img(src='hn_logo.png',
                                   style="margin-top: -10px;
                               padding-right:10px;
                               padding-bottom:10px",
                                   height = 60), "Host & Nested Gene Pair Viewer"),
                   windowTitle="Host & Nested Gene Pair Viewer",
                   
                   # Tab for About Page and additional information to do with the app #    
                   
                   tabPanel("About",
                            fluidRow(
                              column(12,
                                     includeMarkdown("rmds/about.Rmd")),
                            ), 
                            fluidRow(column(4,plotOutput("pipeline"))),
                   ), 
                   
                   
                   
                   
           # Tab for Human #
           
           tabPanel("Human (hg19) Pairs",
                    fluidPage(
                      theme = shinytheme("sandstone"),
                      fluidRow(
                        HTML(
                          paste(
                            h3("Hg19 - Host & Nested Gene Pairs"),'<br/>',
                            p("Table containing list of all pairs of Host & Nested Genes, each respective pair is associated with a pair_id.
                              Click on the row to pull out a specific pair_id which will be used for subsequent plotting and visualisation within the app."),
                            p("If you have a particular or favourite gene, then you can utilise the search function on the top right of the table.")
                          )
                        )
                      ),
                      fluidRow(DT::dataTableOutput("table")),
                      fluidRow(verbatimTextOutput("selection_text")),
                      fluidRow(
                        HTML(
                          paste(
                            '<br/>',
                            p("The two options below will display expression of a selected Host & Nested gene pair across multiple tissues. Normalised counts were calculate from DESeq2."),
                            p("Showing either, the expression of the Host & Nested gene across tissues and, if they are correlated (spearman). More information regarding what ENCODE RNA-seqencing datasets were used can be found in the About tab."),
                            p("A row in the above table must be selected / highlighted for the plotting and visualisation to work. ")
                          )
                        )
                      ),
                      fluidRow(
                        column(4,actionButton("plot_button2", "View expression across tissues")),
                        column(8,actionButton("plot_button", "View spearman correlation scatter plot"))
                      ), 
                      fluidRow(
                        column(4,  plotOutput("expression2")),
                        column(8, plotOutput("expression"))
                      ),
                      fluidRow(
                        HTML(
                          paste(
                            '<br/>',
                            p("Click the button below to generate a customised UCSC link to view the corresponding Host & Nested gene pair. More information regarding what is shown on this custom session can be found at the about page.")
                          )
                        ), 
                        actionButton("url_button", "Generate UCSC Link"),
                        uiOutput("tab"), 
                        HTML(
                          paste(
                            '<br/>')
                        ),
                      ),
                    )), 
           
           
           
           # Tab for Mouse # 
           
           tabPanel("Mouse (mm10) Pairs",
                    fluidPage(
                      theme = shinytheme("sandstone"),
                      fluidRow(
                        HTML(
                          paste(
                            h3("mm10 - Host & Nested Gene Pairs"),'<br/>',
                            p("Table containing list of all pairs of Host & Nested Genes, each respective pair is associated with a pair_id.
                    Click on the row to pull out a specific pair_id which will be used for subsequent plotting and visualisation within the app."),
                            p("If you have a particular or favourite gene, then you can utilise the search function on the top right of the table.")
                          )
                        )
                      ),
                      fluidRow(DT::dataTableOutput("m_table")),
                      fluidRow(verbatimTextOutput("m_selection_text")),
                      fluidRow(
                        HTML(
                          paste(
                            '<br/>',
                            p("The two options below will display expression of a selected Host & Nested gene pair across multiple tissues. Normalised counts were calculate from DESeq2."),
                            p("Showing either, the expression of the Host & Nested gene across tissues and, if they are correlated (spearman). More information regarding what ENCODE RNA-seqencing datasets were used can be found in the About tab."),
                            p("A row in the above table must be selected / highlighted for the plotting and visualisation to work. ")
                          )
                        )
                      ),
                      fluidRow(
                        column(4,actionButton("m_plot_button2", "View expression across tissues")),
                        column(8,actionButton("m_plot_button", "View spearman correlation scatter plot"))
                      ), 
                      fluidRow(
                        column(4,  plotOutput("m_expression2")),
                        column(8, plotOutput("m_expression"))
                      ),
                      fluidRow(
                        HTML(
                          paste(
                            '<br/>',
                            p("Click the button below to generate a customised UCSC link to view the corresponding Host & Nested gene pair. More information regarding what is shown on this custom session can be found at the about page.")
                          )
                        ), 
                        actionButton("m_url_button", "Generate UCSC Link"),
                        uiOutput("m_tab"), 
                        HTML(
                          paste(
                            '<br/>')
                        ),
                      ),
                    )), 
           
           # Tab for Testis Data  # 
           
           tabPanel("Spermatogenesis scRNA-seq",
                    fluidPage(
                      theme = shinytheme("sandstone"),
                      titlePanel("Host / Nested Gene cross talk throughout spermatogenesis"),   
                      sidebarLayout(
                        sidebarPanel(
                          p("Select from the drop down list a testis specific host/nested gene pair to visulise co-expression data."),
                          selectInput("variable", "Select Host:Nested Gene pair", unique(proportion3$H_N)),
                          p("The line plot shows host nested expression proportions across spermatogenesis, the scatter plot shows log(normalised counts) expression of single cells at each of those stages"),
                          
                        ), 
                        mainPanel(
                          plotOutput(outputId = "plot",  width = "100%"),
                          plotOutput(outputId = "plot2", height=650, width = "100%")
                        )
                      )
                    )), 
           
           
           
           tabPanel("Download Data", 
                    fluidRow(column(4,
                                    HTML(
                                      paste(
                                        '<br/>',
                                        '<br/>',
                                        '<br/>',
                                        h3("   ENCODE RNA-Sequencing Metadata"), 
                                        p("   Click to download the metadata from the RNA-sequencing analysis of ENCODE datasets"),
                                        p("")
                                      )),
                                    downloadButton("downloadHumanMetaData", "Human ENCODE Metadata"), 
                                    downloadButton("downloadMouseMetaData", "Mouse ENCODE Metadata")
                    )
                    ), 
                    fluidRow(column(4,
                                    HTML(
                                      paste(
                                        '<br/>',
                                        '<br/>',
                                        '<br/>',
                                        h3("   The Host / Nested Gene Data"), 
                                        p("   Click to download the dataset generated in Montibus.et al containing all Human and Mouse Host/Nested Genes"),
                                        p("")
                                      )),
                                    downloadButton("downloadHumanHnGenes", "Hg19 Host Nested Genes"), 
                                    downloadButton("downloadMouseHnGenes", "Mm10 Host Nested Genes")
                    )
                    ), 
           )
)
)


