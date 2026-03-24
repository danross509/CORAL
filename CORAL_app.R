#
# Shiny application to filter and generate interactive graphics for DESeq2 results 
#
# Requires:   Count matrix as csv or tsv file
#             Metadata as csv or tsv file, corresponding to the samples in the count matrix
#             EITHER folder containing DESeq2 analysis results (including padj and log2foldchange)
#               OR a gene list file containing a single column of gene names as a csv or tsv
#
# Generates:  Filtered DESeq2 results based on input and specified filters
#             Unique gene lists based on specified filters
#             Heatmap of all gene counts, and cluster specific heatmaps
#             GO term plots, if given custom GO terms
#
# Created by: David Ross
#             For the Marine Climate Change Unit at the Okinawa Institute of Science and Technology
#             March 2026
#
# Please use CORAL_requirements.Rmd to ensure required packages are installed before running
#
options(shiny.maxRequestSize = 20 * 1024^2)

library("ggplot2") 
library("ggfortify")
library("gprofiler2")
library("ggrepel")
library("RColorBrewer")
library("pheatmap")
library("EnhancedVolcano")
library("gplots")
library("DESeq2")
library("shinyjs")
library("devtools")
library("colourpicker")
library("data.table")
library("DT")
library("clusterProfiler")
library("enrichplot")
library("DOSE")
library("rrcov")
library("miniUI")
library("shinyFiles")
library("stringr")
library("cluster")
library("plyr")
library("dplyr")

base_plot_colours <- c("#481c6e","#2ab07f","#eae51a")

shinyApp(
  ui = tagList(
    useShinyjs(),
    navbarPage("DESeq2 Downstream Analysis",
               
      ########################################
      #    Data input - UI                   # 
      ########################################               
               
      tabPanel("Data Input",
               h4("Input the requested files"),
               
               sidebarPanel(
                 fileInput("countFile", "TPM gene count file", accept = c(".tsv", ".csv")),
                 div(style = "margin-top: -30px"),
                 p("Set gene count read.table parameters:"),
                 checkboxInput("countHeader", "Count file header", TRUE),
                 checkboxInput("countRownames", "Use column 1 as row names", TRUE),
                 checkboxInput("countSalmon", "Count data generated with Salmon (tx + gene_id)", FALSE),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 fileInput("metaData", "Corresponding metadata file", accept = c(".tsv", ".csv")),
                 div(style = "margin-top: -30px"),
                 p("Set metadata read.table parameters:"),
                 checkboxInput("metaDataHeader", "Metadata file header", TRUE),
                 checkboxInput("metaDataRownames", "Use column 1 as row names", TRUE),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 fileInput("uniqueGeneListInput", "Unique gene list", accept = c(".tsv", ".csv")),
                 div(style = "margin-top: -30px"),
                 checkboxInput("uniqueGeneListHeader", "Unique gene list file header", TRUE),
                 
                 p(strong("OR")),
                 
                 p("Choose a directory:"),
                 shinyDirButton(
                   id = "deseq2FilesDir", 
                   label = "Original DESeq2 files",
                   title = "Choose a directory:"),
                 p("Alternate directories (C:, G:, etc) can be chosen from the dropdown menu"),
                 checkboxInput("deseq2FilesHeader", "DESeq2 file headers", TRUE),
                 checkboxInput("deseq2FilesRownames", "Use column 1 as row names", TRUE),
                 uiOutput("deseq2FilesSelect"),
                 textInput("deseq2padjCol", "DESeq2 p-adj column name:", value="padj"),
                 textInput("deseq2LfcCol", "DESeq2 L2FC column name:", value="log2FoldChange")
                 ),
               
               mainPanel(
                 tabsetPanel(
                  tabPanel("Gene counts", dataTableOutput("countsTable")),
                  tabPanel("Sample metadata", dataTableOutput("metaDataTable")),
                  tabPanel("Unique DEGs", dataTableOutput("uniqueGeneListInputTable")),
                  tabPanel("Original DESeq2 files", verbatimTextOutput("deseq2FileList"))
                 )
               )
                
              ),
      ########################################
      #     Data processing - UI             #
      ########################################
      
      tabPanel("Data processing",
               h4("Choose parameters to filter data"),
               sidebarPanel(
                 p("Select a parent directory to begin processing data"),
                 shinyDirButton(
                    id = "deseq2FilesDirParent", 
                    label = "Parent directory for filtered files", 
                    title = "Choose a directory"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                  
                 numericInput("padj_threshold",
                           "Set adjusted p-value maximum",
                           min = 0,
                           max = 1,
                           value = 0.05,
                           step = 0.005),
                 numericInput("lfc_threshold",
                           "Set log2 fold change minimum",
                           min = 0,
                           max = NA,
                           value = 0.58,
                           step = 0.02),
                 uiOutput("uniqueGeneListDir"),
                 textInput("uniqueGeneListFilePrefix","Prefix for unique DEG list:",value="full_unique_DEGlist"),
                 textInput("deseq2FilesDirFiltPrefix","Prefix for filtered files directory:",value="filtered"),
                 checkboxInput("overwriteUniqueGeneListFile", "Overwrite identical unique DEG files", FALSE),
                 checkboxInput("overwriteDeseq2FilesFilt", "Overwrite identical filtered DESeq2 files", FALSE),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 actionButton("saveFilteredData",
                            "Save filtered data")
               ),
               
               mainPanel(
                 textOutput("uniqueGeneListMessage"),
                 verbatimTextOutput("uniqueGeneListPath"),
                 textOutput("uniqueGeneListWarningOverwriteFiles"),
                 verbatimTextOutput("uniqueGeneListExistingFiles"),
                 textOutput("filteredSaveMessage"),
                 verbatimTextOutput("deseq2FilesDirFiltPath"),
                 textOutput("deseq2WarningOverwriteFiles"),
                 verbatimTextOutput("deseq2ExistingFiltFiles"),
                 textOutput("deseq2FilesListFiltMessage"),
                 verbatimTextOutput("deseq2FilesListFilt")
               )
               
               
        
      ),
      ########################################
      #     Heatmap - UI                     #
      ########################################      
      
      tabPanel("Heatmap",         
             h4("Heatmap (Please wait for data to process)"),     
             sidebarPanel(
               strong("Data input options:"),
               checkboxInput("hm_useGroupMean", "Use per-treatment mean counts", TRUE),
               checkboxInput("hm_useGeneListAuto", "Prioritize auto-generated gene list from original DESeq2 data", TRUE),
               selectizeInput("hm_scaling",
                              "Select data scaling metric",
                              c("absolute TPM", "log2(TPM+1)", "z-score", "z-score(log2(TPM+1))"),
                              selected="log2(TPM+1)"),
               tags$hr(style="background-color: black; height: 1px; border: 0"),
               
               strong("Clustering options:"),
               numericInput("hm_geneClusterK", 
                            "Choose # of clusters to generate",
                            min=2,
                            max=10,
                            value=3),
               actionButton("hm_showGeneClusters", "Show gene clusters"),
               div(style = "margin-top: +10px"),
               textInput("hm_clusterFilePrefix","Prefix for per-cluster DEG files:",value="cluster_DEGs"),
               actionButton("hm_saveGeneClusterFiles", "Save per-cluster DEG files"),
               div(style = "margin-top: +10px"),
               checkboxInput("saveMainClusterTxt", "Also save cluster files as .txt", FALSE),
               checkboxInput("overwriteMainClusterFiles", "Overwrite identical cluster files", FALSE),
               selectizeInput("hm_showDendrogram", "Clustering:",
                              choices=c("both" = "both",
                                        "on samples" = "column",
                                        "on genes" = "row",
                                        "none" = "none"),
                              selected="both",
                              multiple=FALSE),
               selectizeInput("hm_cluster_dist",
                              "Clustering distance measure:",
                              choices=c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                              selected="correlation",
                              multiple=FALSE),
               uiOutput("hm_cutreeCols"),
               selectizeInput("hm_cutreeRows", 
                              "Split genes into X clusters:", 
                              choices=c("None", sequence(12, from=2L)), 
                              selected="None", 
                              multiple=FALSE),
               tags$hr(style="background-color: black; height: 1px; border: 0"),
               
               strong("Esthetic options:"),
               uiOutput("hm_title_ui"),
               selectizeInput("hm_color",
                              "Select ColorBrewer palette:",
                              c("Oranges", "Blues", "Reds", "BuGn", "Greens", "Purples", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"),
                              selected="BuGn"
               ),
               checkboxInput("hm_color_rev", "Reverse color palette:", FALSE),
               sliderInput("hm_color_number",
                           "Select number of distinct colors:",
                           min = 3,
                           max = 12,
                           value = 9,
                           step = 1),
               sliderInput("hm_base_fontsize",
                           "Adjust plot base fontsize",
                           min = 0.1,
                           max = 20,
                           value = 10,
                           step = 0.1),
               checkboxInput("hm_gene_names", "Show gene names:", FALSE),
               sliderInput("hm_gene_fontsize",
                           "Gene name font size:",
                           min = 0.1,
                           max = 5,
                           value = 0.8,
                           step = 0.1),
               checkboxInput("hm_sample_names", "Show sample names:", TRUE),
               selectizeInput("hm_sampleName_angle",
                              "Choose sample name display angle:",
                              c("270", "0", "45", "90", "315"),
                              selected="270"),
               sliderInput("hm_sample_fontsize",
                           "Sample name font size:",
                           min = 0.1,
                           max = 20,
                           value = 10,
                           step = 0.5),
               checkboxInput("hm_showLegend", "Show legend:", TRUE),
               textInput("hm_legenTitle", "Colour key title", value=""),
               colourInput("hm_sepcolor", "Select separation colour", "white"),
               tags$hr(style="background-color: black; height: 1px; border: 0"),
               
               strong("Image download options:"),
               uiOutput("hm_filename"),
               selectInput("hm_file_format", "File Format", 
                           choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                           selected="png"),
               numericInput("hm_plot_width", "Width (inches)", value = 12, min = 1, max = 50),
               numericInput("hm_plot_height", "Height (inches)", value = 18, min = 1, max = 50),
               numericInput("hm_plot_dpi", "Resolution (DPI)", value = 300, min = 72, max = 600)
             ),
             
             mainPanel(plotOutput("HeatmapPlot", width = "80%", height = "800px"),
                       downloadButton("download_heatmap", "Download heatmap"),
                       div(style = "margin-top: +10px"),
                       p("Cluster files will be saved under the following directory:"),
                       verbatimTextOutput("geneClusterFilesPath"),
                       textOutput("clusterWarningOverwriteFiles"),
                       verbatimTextOutput("clusterMainExistingFiles"),
                       actionButton("hm_generateClusterSil", "Show cluster silhouettes"),
                       plotOutput("hm_clusterSilPlot", width = "80%")
                       ),
            ),
      
      ########################################
      #     Sub-cluster heatmap - UI    #
      ########################################    
      tabPanel("Sub-cluster heatmap",         
               h4("Cluster-specific heatmap"),     
               sidebarPanel(
                 uiOutput("hmSUB_clusterChoice"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Clustering options:"),
                 numericInput("hmSUB_geneClusterK", 
                              "Choose # of clusters to generate",
                              min=2,
                              max=10,
                              value=2),
                 actionButton("hmSUB_showGeneClusters", "Show gene clusters"),
                 textInput("hmSUB_clusterFilePrefix","Prefix for sub-cluster filenames:",value="sub-DEGs"),
                 actionButton("hmSUB_saveGeneClusterFiles", "Save sub-cluster DEG files"),
                 div(style = "margin-top: +10px"),
                 checkboxInput("saveSUBClusterTxt", "Also save sub-cluster files as .txt", FALSE),
                 checkboxInput("overwriteSUBClusterFiles", "Overwrite identical sub-cluster files", FALSE),
                 selectizeInput("hmSUB_showDendrogram", "Clustering:",
                                choices=c("both" = "both",
                                          "on samples" = "column",
                                          "on genes" = "row",
                                          "none" = "none"),
                                selected="row",
                                multiple=FALSE),
                 uiOutput("hmSUB_cutreeCols"),
                 selectizeInput("hmSUB_cutreeRows", 
                                "Split genes into X clusters:", 
                                choices=c("None", sequence(12, from=2L)), 
                                selected="None", 
                                multiple=FALSE),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 uiOutput("hmSUB_title_ui"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Image download options:"),
                 uiOutput("hmSUB_filename"),
                 selectInput("hmSUB_file_format", "File Format", 
                             choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                             selected="png"),
                 numericInput("hmSUB_plot_width", "Width (inches)", value = 12, min = 1, max = 50),
                 numericInput("hmSUB_plot_height", "Height (inches)", value = 18, min = 1, max = 50),
                 numericInput("hmSUB_plot_dpi", "Resolution (DPI)", value = 300, min = 72, max = 600)
               ),
               
               mainPanel(plotOutput("SUBHeatmapPlot", width = "80%", height = "800px"),
                         downloadButton("download_SUBheatmap", "Download sub-cluster heatmap"),
                         p("Sub-cluster files will be saved under the following directory:"),
                         verbatimTextOutput("geneSUBClusterFilesPath"),
                         textOutput("SUBclusterWarningOverwriteFiles"),
                         verbatimTextOutput("SUBclusterMainExistingFiles"),
                         actionButton("hmSUB_generateClusterSil", "Show cluster silhouettes"),
                         plotOutput("hmSUB_clusterSilPlot", width = "80%")
                         )
            ),

      ########################################
      #     Sub-sub-cluster heatmap - UI #
      ########################################
      tabPanel("Sub2-cluster heatmap",         
               h4("Sub2-cluster heatmap"),     
               sidebarPanel(
                 uiOutput("hmSUB2_clusterChoice"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Clustering options:"),
                 numericInput("hmSUB2_geneClusterK", 
                              "Choose # of clusters to generate",
                              min=2,
                              max=10,
                              value=2),
                 actionButton("hmSUB2_showGeneClusters", "Show gene clusters"),
                 textInput("hmSUB2_clusterFilePrefix","Prefix for sub2-cluster filenames:",value="sub2-DEGs"),
                 actionButton("hmSUB2_saveGeneClusterFiles", "Save sub2-cluster DEG files"),
                 div(style = "margin-top: +10px"),
                 checkboxInput("saveSUB2ClusterTxt", "Also save sub-cluster files as .txt", FALSE),
                 checkboxInput("overwriteSUB2ClusterFiles", "Overwrite identical sub-cluster files", FALSE),
                 selectizeInput("hmSUB2_showDendrogram", "Clustering:",
                                choices=c("both" = "both",
                                          "on samples" = "column",
                                          "on genes" = "row",
                                          "none" = "none"),
                                selected="row",
                                multiple=FALSE),
                 uiOutput("hmSUB2_cutreeCols"),
                 selectizeInput("hmSUB2_cutreeRows", 
                                "Split genes into X clusters:", 
                                choices=c("None", sequence(12, from=2L)), 
                                selected="None", 
                                multiple=FALSE),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 uiOutput("hmSUB2_title_ui"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Image download options:"),
                 uiOutput("hmSUB2_filename"),
                 selectInput("hmSUB2_file_format", "File Format", 
                             choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                             selected="png"),
                 numericInput("hmSUB2_plot_width", "Width (inches)", value = 12, min = 1, max = 50),
                 numericInput("hmSUB2_plot_height", "Height (inches)", value = 18, min = 1, max = 50),
                 numericInput("hmSUB2_plot_dpi", "Resolution (DPI)", value = 300, min = 72, max = 600)
               ),
               
               mainPanel(plotOutput("SUB2HeatmapPlot", width = "80%", height = "800px"),
                         downloadButton("download_SUB2heatmap", "Download heatmap"),
                         p("Sub2-cluster files will be saved under the following directory:"),
                         verbatimTextOutput("geneSUB2ClusterFilesPath"),
                         textOutput("SUB2clusterWarningOverwriteFiles"),
                         verbatimTextOutput("SUB2clusterMainExistingFiles"),
                         actionButton("hmSUB2_generateClusterSil", "Show cluster silhouettes"),
                         plotOutput("hmSUB2_clusterSilPlot", width = "80%")
                        )
               ),
      
      ########################################
      #     Cluster GO-terms - UI            #
      ########################################
      tabPanel("Cluster GO-terms",         
               h4("Cluster GO-terms"),     
               sidebarPanel(
                 fileInput("go2geneConversion", "GO:ID <-> gene conversion file", accept = c(".tsv", ".csv", ".txt")),
                 div(style = "margin-top: -30px"),
                 p("Set GO:ID <-> gene read.table parameters:"),
                 checkboxInput("go2geneHeader", "GO:ID <-> gene file header", FALSE),
                 checkboxInput("go2geneRownames", "Use column 1 as row names", FALSE),
                 uiOutput("go2geneTxtExt"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 fileInput("go_termConversion", "GO:ID <-> term name conversion file", accept = c(".tsv", ".csv", ".txt")),
                 div(style = "margin-top: -30px"),
                 p("Set GO:ID <-> term read.table parameters:"),
                 checkboxInput("go2termHeader", "GO:ID <-> term file header", FALSE),
                 checkboxInput("go2termRownames", "Use column 1 as row names", FALSE),
                 uiOutput("go2termTxtExt"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 selectizeInput("go_hmChoice", "Choose heatmap level:",
                                choices=c("Main",
                                          "Sub",
                                          "Sub2"),
                                selected="Sub2",
                                multiple=FALSE),
                 uiOutput("go_clusterChoice"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Set ClusterProfiler::enricher parameters:"),
                 p("(padj and q value limits will be set in the bubble plot tab)"),
                 selectizeInput("ego_pAdjustMethod",
                                "Enricher pAdjustMethod",
                                choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                selected="BH",
                                multiple=FALSE),
                 numericInput("ego_minGSSize",
                              "Enricher minGSSize",
                              min = 0,
                              max = 1000,
                              value = 5,
                              step = 1),
                 numericInput("ego_maxGSSize",
                              "Enricher maxGSSize",
                              min = 0,
                              max = NA,
                              value = 5000,
                              step = 50),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 p("Modify 'preview' plot output:"),
                 selectizeInput("go_largePlotStyle",
                                "Large plot display style",
                                choices=c("bubble", "bar"),
                                selected="bubble",
                                multiple=FALSE),
                 textInput("go_largePlotTitle", "Large plot title", value="GO enrichment (custom annotation)"),
                 uiOutput("go_previewLimit"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 p("Modify histogram output:"),
                 uiOutput("go_countLimit"),
                 
                 strong("'Preview' image download options:"),
                 uiOutput("goPreview_filename"),
                 selectInput("goPreview_file_format", "File Format", 
                             choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                             selected="png"),
                 numericInput("goPreview_plot_width", "Width (inches)", value = 12, min = 1, max = 50),
                 numericInput("goPreview_plot_height", "Height (inches)", value = 18, min = 1, max = 50),
                 numericInput("goPreview_plot_dpi", "Resolution (DPI)", value = 300, min = 72, max = 600)
               ),

               mainPanel(plotOutput("go_previewPlot", width = "80%", height = "800px"),
                         downloadButton("download_goPreview", "Download 'preview' plot"),
                         plotOutput("go_countHistogram", width = "80%", height = "800px")
                         )
              ),
      
      ########################################
      #     Cluster bubble plot - UI         #
      ########################################
      tabPanel("Cluster bubble plot",
               h4("Cluster bubble plot"),
               sidebarPanel(
                 strong("Modify bubble plot data selection parameters:"),
                 numericInput("ego_pvalueCutoff",
                              "Enricher pvalueCutoff",
                              min = 0,
                              max = 1,
                              value = 0.05,
                              step = 0.005),
                 numericInput("ego_qvalueCutoff",
                              "Enricher qvalueCutoff",
                              min = 0,
                              max = 1,
                              value = 1,
                              step = 0.01),
                 selectizeInput("ego_dfSort",
                                "Arrange GO db values (for cutoff) by:",
                                choices=c("RichFactor",
                                          "FoldEnrichment",
                                          "zScore",
                                          "p.adjust",
                                          "qvalue",
                                          "Count",
                                          "GeneRatio_num",
                                          "neglogFDR"),
                                selected="neglogFDR",
                                multiple=FALSE),
                 checkboxInput("ego_dfSortDesc", "Sort descending:", TRUE),
                 uiOutput("go_bubbleLimit"),
                 checkboxInput("go_bubbleChooseTerms", "Manually select GO terms to display?", FALSE),
                 uiOutput("go_bubbleTermChoice"),
                 tags$hr(style="background-color: black; height: 1px; border: 0"),
                 
                 strong("Modify bubble plot display parameters:"),
                 p("Bubble parameters"),
                 selectizeInput("ego_plotSort",
                                "Arrange GO bubble plot display by:",
                                choices=c("RichFactor",
                                          "FoldEnrichment",
                                          "zScore",
                                          "p.adjust",
                                          "qvalue",
                                          "Count",
                                          "GeneRatio_num",
                                          "neglogFDR"),
                                selected="Count",
                                multiple=FALSE),
                 checkboxInput("ego_plotSortDesc", "Sort descending:", TRUE),
                 colourInput("go_bubbleColorHigh", "Select 'high' colour", base_plot_colours[3]),
                 colourInput("go_bubbleColorMid", "Select 'mid' colour", base_plot_colours[2]),
                 colourInput("go_bubbleColorLow", "Select 'low' colour", base_plot_colours[1]),
                 checkboxInput("ego_revPlotColours", "Reverse colours", FALSE),
                 checkboxInput("go_bubbleDotOutline", "Outline bubbles", TRUE),
                 numericInput("ego_geomOutlineSize",
                              "Point outline thickness",
                              min = 0,
                              max = 10,
                              value = 1,
                              step = 0.2),
                 numericInput("ego_geomAlpha",
                              "Point color opacity",
                              min = 0,
                              max = 1,
                              value = 0.85,
                              step = 0.05),
                 
                 p("General settings"),
                 uiOutput("go_bubbleTitle"),
                 numericInput("ego_titleSize",
                              "Plot title size",
                              min = 0,
                              max = 40,
                              value = 12,
                              step = 1),
                 selectizeInput("ego_bubblePlotTheme",
                                "Bubble plot theme",
                                choices=c("bw",
                                          "linedraw",
                                          "light",
                                          "dark",
                                          "minimal",
                                          "classic",
                                          "gray"),
                                selected="bw",
                                multiple=FALSE),
                 numericInput("ego_themeBaseSize",
                              "Theme base line size",
                              min = 2,
                              max = 40,
                              value = 12,
                              step = 2),
                 
                 p("Y-axis label parameters"),
                 numericInput("ego_yAxisStrWrap",
                              "Y-axis label wrap length",
                              min = 5,
                              max = 100,
                              value = 20,
                              step = 5),
                 checkboxInput("go_yAxisShowTitle", "Show y-axis title", FALSE),
                 textInput("go_yAxisTitle", "y-axis title", value=""),
                 numericInput("ego_yAxisTextSize",
                              "Y-axis text size",
                              min = 0,
                              max = 40,
                              value = 9,
                              step = 1),
                 checkboxInput("go_yAxisTicks", "Show y-axis ticks", FALSE),
                 checkboxInput("go_yAxisGridlines", "Show y-axis (horizontal) gridlines", TRUE),
                 
                 p("X-axis label parameters"),
                 checkboxInput("go_xAxisShowTitle", "Show x-axis title", TRUE),
                 textInput("go_xAxisTitle", "X-axis title", value="Gene Ratio"),
                 numericInput("ego_xAxisTextSize",
                              "X-axis text size",
                              min = 0,
                              max = 40,
                              value = 9,
                              step = 1),
                 checkboxInput("go_xAxisTicks", "Show x-axis ticks", TRUE),
                 checkboxInput("go_xAxisGridlines", "Show x-axis (vertical) gridlines", TRUE),
                 
                 p("Legend parameters"),
                 checkboxInput("ego_legendBox", "Frame legend", TRUE),
                 selectizeInput("ego_bubbleLegendPosition",
                                "Legend position",
                                choices=c("none", 
                                          "left", 
                                          "right", 
                                          "bottom", 
                                          "top", 
                                          "inside"),
                                selected="right",
                                multiple=FALSE),
                 numericInput("ego_legendSpacing",
                              "Legend spacing (outside plot)",
                              min = 0,
                              max = 40,
                              value = 10,
                              step = 1),
                 selectizeInput("ego_bubbleLegendPositionInside",
                                "Legend alignment (inside)",
                                choices=c("bottom left", 
                                          "bottom right", 
                                          "top left", 
                                          "top right"),
                                selected="bottom right",
                                multiple=FALSE),
                 selectizeInput("ego_bubbleLegendStacking",
                                "Stack legends:",
                                choices=c("vertical", 
                                          "horizontal"),
                                selected="vertical",
                                multiple=FALSE),
                 numericInput("ego_legendMargin",
                              "Legend margin",
                              min = 0,
                              max = 40,
                              value = 6,
                              step = 1),
                 numericInput("ego_legendTitleSize",
                              "Legend title size",
                              min = 0,
                              max = 40,
                              value = 10,
                              step = 1),
                 numericInput("ego_legendTextSize",
                              "Legend text size",
                              min = 0,
                              max = 40,
                              value = 9,
                              step = 1),
                 
                 strong("Image download options:"),
                 uiOutput("goEGO_filename"),
                 selectInput("goEGO_file_format", "File Format", 
                             choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg", "TIFF" = "tiff"),
                             selected="png"),
                 numericInput("goEGO_plot_width", "Width (inches)", value = 12, min = 1, max = 50),
                 numericInput("goEGO_plot_height", "Height (inches)", value = 18, min = 1, max = 50),
                 numericInput("goEGO_plot_dpi", "Resolution (DPI)", value = 300, min = 72, max = 600)
               ),

               mainPanel(plotOutput("go_EGO_plot", width = "80%", height = "800px"),
                         downloadButton("download_goEGO", "Download bubble plot"))
            )
      )
    ),
   
    server = function(input, output, session) {
    
    options(shiny.maxRequestSize=100*1024^2)
    
    ########################################
    #    Data input - Server               #
    ########################################
    
    # Input raw TPM data file
    countData <- reactive({
      req(input$countFile)
      countFile <- input$countFile
      countExt <- tools::file_ext(countFile$datapath)
      
      validate(need(countExt %in% c("tsv", "csv"), "Please upload either a csv or tsv file"))
      if (countExt == "tsv") {
        countSep <- "\t"
      } else if (countExt == "csv") {
        countSep <- ","
      }
      
      if (input$countRownames) {
        countRownames <- 1
      } else {
        countRownames <- NULL
      }
      
      read.table(
        countFile$datapath, 
        header=input$countHeader, 
        row.names = countRownames, 
        sep=countSep, 
        fill=FALSE)
    })  
      
    output$countsTable <- renderDataTable({
      countData()
    })
  
    
    # Input sample meta data file
    metaData <- reactive({
      req(input$metaData)
      metaData <- input$metaData
      metaDataExt <- tools::file_ext(metaData$datapath)
      
      validate(need(metaDataExt %in% c("tsv", "csv"), "Please upload either a csv or tsv file"))
      if (metaDataExt == "tsv") {
        metaSep <- "\t"
      } else if (metaDataExt == "csv") {
        metaSep <- ","
      }
      
      if (input$metaDataRownames) {
        metaDataRownames <- 1
      } else {
        metaDataRownames <- NULL
      }
      
      read.table(
        metaData$datapath, 
        header=input$metaDataHeader, 
        row.names = metaDataRownames,
        colClasses = "character",
        sep=metaSep, 
        fill=FALSE)
    })  
    
    output$metaDataTable <- renderDataTable({
      metaData()
    })
    
    uniqueGeneListInputTable <- reactive ({
      req(input$uniqueGeneListInput)
      
      uniqueDEGs <- input$uniqueGeneListInput
      uniqueDEGsExt <- tools::file_ext(uniqueDEGs$datapath)
      
      validate(need(uniqueDEGsExt %in% c("tsv", "csv"), "Please upload either a csv or tsv file"))
      if (uniqueDEGsExt == "tsv") {
        DEGSep <- "\t"
      } else if (uniqueDEGsExt == "csv") {
        DEGSep <- ","
      }
      
      read.table(
        uniqueDEGs$datapath, 
        header=input$uniqueGeneListHeader,
        row.names=NULL,
        sep=DEGSep, 
        fill=FALSE)
    })
    
    output$uniqueGeneListInputTable <- renderDataTable({
      uniqueGeneListInputTable()
    })
    
    
    # Input DESeq2 file directory
    volumes <- c(Home = fs::path_home(), getVolumes()())
    shinyDirChoose(input, 'deseq2FilesDir', roots = volumes, session = session)
    
    deseq2FilesDir <- reactive({
      req(input$deseq2FilesDir)
      
      parseDirPath(volumes, input$deseq2FilesDir)
    })
    
    output$deseq2FilesDir <- renderPrint(deseq2FilesDir())
    
    # Choose files, if necessary
    
    output$deseq2FilesSelect <- renderUI({
      req(deseq2FilesDir())
      
      path <- deseq2FilesDir()
      
      choices <- list.files(path, full.names = FALSE)
      
      selectizeInput("deseq2FilesSelect",
                     "Select DESeq2 files to use:",
                     c(choices),
                     selected=choices,
                     multiple=TRUE)
    })
    
    # Path to chosen files
    output$deseq2FileList <- renderPrint({
      req(deseq2FilesDir())
      
      path <- deseq2FilesDir()
      
      deseq2FileList <- list.files(path, full.names = TRUE)
      
      finalFileList <- deseq2FileList[which(basename(deseq2FileList) %in% input$deseq2FilesSelect)]
      
      if (length(finalFileList) == 0) {
        cat("Directory is empty")
      } else {
        cat(paste(finalFileList, collapse = "\n"))
      }
    })

    ########################################
    #     Data processing - Server         #
    ########################################
    
    # Display path for new file directory
    output$filteredSaveMessage <- renderText("Filtered DESeq2 result files will be saved in the directory:")
    output$uniqueGeneListMessage <- renderText("Unique DEG lists will be saved in the directory:")
    output$deseq2FilesListFiltMessage <- renderText("Filtered DESeq2 files to be saved:")
    
    # Choose directory for filtered DESeq2 files
    parent_DESeq2FilesDir <- reactive({
      req(deseq2FilesDir())
      
      parent_dir <- dirname(deseq2FilesDir())
      
      c(Home = parent_dir, getVolumes()())
    })
    
    observe({
      shinyDirChoose(input, 
                     'deseq2FilesDirParent', 
                     roots = parent_DESeq2FilesDir(), 
                     session = session)
    })
    
    deseq2FilesDirParent <- reactive({
      req(input$deseq2FilesDirParent)
      
      parseDirPath(parent_DESeq2FilesDir(), input$deseq2FilesDirParent)
    })
    
    # Create directory name based on prefix + selected thresholds
    output$deseq2FilesDirFiltPath <- renderPrint(paste(deseq2FilesDirParent(), 
                                                   "/", 
                                                   input$deseq2FilesDirFiltPrefix, 
                                                   "_padj", 
                                                   padj_d(), 
                                                   "_log2FC", 
                                                   lfc_d(),
                                                   "/", sep=""))
    
    # Set unique_DEG directory name based on input
    output$uniqueGeneListDir <- renderUI({
      req(input$deseq2FilesSelect)
      
      if (length(input$deseq2FilesSelect) > 1) {
        value <- "combined_unique_DEGlists"
      } else {
        value <- "unique_DEGlists"
      }
      
      textInput("uniqueGeneListDir",
                "Name for unique DEG directory:",
                value=value)
    })
    
    # Create directory for unique gene list
    output$uniqueGeneListPath <- renderPrint(paste(deseq2FilesDirParent(), 
                                                   "/", input$uniqueGeneListDir, "/", 
                                                   sep=""))
    
    # If gene list folder exists and contains files, display warning
    uniqueGeneListFile <- reactive({
      req(deseq2FilesDirParent())
      
      path <- paste(deseq2FilesDirParent(), 
                    "/", input$uniqueGeneListDir, "/",
                    input$uniqueGeneListFilePrefix,
                    "_padj", 
                    padj_d(), 
                    "_log2FC", 
                    lfc_d(),
                    ".csv",
                    sep="")
      
      return(path)
    })
      
    output$uniqueGeneListWarningOverwriteFiles <- renderText({
      req(deseq2FilesDirParent())
      update_saveFilteredData()
      
      uniqueDEGFile <- uniqueGeneListFile()
      
      if (!file.exists(uniqueDEGFile)) {
        return("Unique DEG file not yet generated with given parameters")
      } else if (input$overwriteUniqueGeneListFile) {
        return("WARNING: Unique DEG directory currently contains the following file. Identical filenames will be overwritten")
      } else {
        return("WARNING: Unique DEG directory currently contains the following file. Identical filenames will be ignored")
      }
    })
    
    # Display current unique gene list files
    output$uniqueGeneListExistingFiles <- renderPrint({
      req(deseq2FilesDirParent())
      update_saveFilteredData()
      
      uniqueDEGFile <- uniqueGeneListFile()
      
      if (!file.exists(uniqueDEGFile)) {
        return(paste("File", uniqueDEGFile, "does not yet exist", sep = " "))
      } else {
        return(uniqueDEGFile)
      }
    })
    
    # If filtered DESeq2 folder exists and contains files, display warning
    deseq2FileListFilt <- reactive({
      req(deseq2FilesDirParent())
      update_saveFilteredData()
      
      path <- paste(deseq2FilesDirParent(), 
                    "/", 
                    input$deseq2FilesDirFiltPrefix, 
                    "_padj", 
                    padj_d(), 
                    "_log2FC", 
                    lfc_d(),
                    "/", sep="")
      
      list.files(path, full.names = FALSE)
    })
    
    output$deseq2WarningOverwriteFiles <- renderText({
      req(deseq2FilesDirParent())
      
      existingFileList <- deseq2FileListFilt()
      
      if (length(existingFileList) == 0) {
        return("Filtered DESeq2 file directory is currently empty")
      } else if (input$overwriteDeseq2FilesFilt) {
        return("WARNING: Filtered DESeq2 file directory currently contains the following files. Identical filenames will be overwritten")
      } else {
        return("WARNING: Filtered DESeq2 file directory currently contains the following files. Identical filenames will be ignored")
      }
    })
    
    # Display current filtered DESeq2 files
    output$deseq2ExistingFiltFiles <- renderPrint({
      req(deseq2FilesDirParent())
      
      existingFileList <- deseq2FileListFilt()
      
      if (length(existingFileList) == 0) {
        return("No files to display")
      } else {
        return(cat(paste(existingFileList, collapse = "\n")))
      }
    })
    
    # Establish DESeq2 filtering input reactive
    padj_d  <- debounce(reactive(input$padj_threshold), 800)
    lfc_d   <- debounce(reactive(input$lfc_threshold),  800)
    
    deseq2InputsReady <- reactive ({
      req(deseq2FilesDir())
      req(padj_d())
      req(lfc_d())
      
      TRUE
    })
    
    # Display list of files to be saved: filtered DESeq2 files AND genelist
    output$deseq2FilesListFilt <- renderPrint({
      req(deseq2InputsReady())
      req(deseq2FilesDirParent())
      
      pathRaw <- deseq2FilesDir()
      pathFilt <- paste(deseq2FilesDirParent(), 
                        "/", 
                        input$deseq2FilesDirFiltPrefix, 
                        "_padj", 
                        padj_d(), 
                        "_log2FC", 
                        lfc_d(),
                        "/", sep="")
      
      deseq2FileList <- list.files(pathRaw, full.names = TRUE)
      finalFileList <- deseq2FileList[which(basename(deseq2FileList) %in% input$deseq2FilesSelect)]
      
      deseq2FilesListFilt <- c()
      for (rawFile in finalFileList) {
        fileExt <- tools::file_ext(rawFile)
        fileName <- tools::file_path_sans_ext(rawFile)
        filtName <- paste(fileName, 
                          "_padj", padj_d(), 
                          "_log2FC", lfc_d(), 
                          ".", fileExt, 
                          sep = "")
        deseq2FilesListFilt <-c(deseq2FilesListFilt, filtName)
      }
      
      
      if (length(deseq2FilesListFilt) == 0) {
        cat("Original directory is empty")
      } else {
        cat(paste(deseq2FilesListFilt, collapse = "\n"))
      }
    })
    
    # Generate unique gene list from input DESeq2 files
    uniqueGeneListAuto <- reactive ({
      req(deseq2InputsReady())
      
      pathRaw <- deseq2FilesDir()
      deseq2FileList <- list.files(pathRaw, full.names = TRUE)
      finalFileList <- deseq2FileList[which(basename(deseq2FileList) %in% input$deseq2FilesSelect)]
      
      if (input$deseq2FilesRownames) {
        deseq2FilesRownames <- 1
      } else {
        deseq2FilesRownames <- NULL
      }
      
      
      allDEGs <- c()
      
      for (rawFile in finalFileList) {
        
        fileExt <- tools::file_ext(rawFile)
        
        if (fileExt == "tsv") {
          fileSep <- "\t"
        } else if (fileExt == "csv") {
          fileSep <- ","
        }
        
        rawFile_df <- read.table(file = rawFile, 
                                 header=input$deseq2FilesHeader, 
                                 row.names = deseq2FilesRownames,
                                 sep=fileSep, 
                                 fill=FALSE)
        
        genesToKeep <- which(!is.na(rawFile_df[, input$deseq2padjCol]) 
                             & rawFile_df[, input$deseq2padjCol]<=padj_d()
                             & abs(as.numeric(rawFile_df[, input$deseq2LfcCol]))>=lfc_d())
        filtFile_df <- rawFile_df[genesToKeep,]
        
        allDEGs <- c(allDEGs, rownames(filtFile_df))
      }
      
      uniqueDEGs <- unique(allDEGs)
      
      return(data.frame(gene_id = uniqueDEGs))
    })
    
    # Reactive value to force reactive updates on file save
    update_saveFilteredData <- reactiveVal(0)
    
    # Generate DESeq2 files and save with unique gene list on button click
    observeEvent(input$saveFilteredData, {
      req(deseq2InputsReady())
      req(deseq2FilesDirParent())
      
      uniqueGeneListDir <- paste(deseq2FilesDirParent(), 
                                 "/", input$uniqueGeneListDir, "/", 
                                 sep="")
      
      pathRaw <- deseq2FilesDir()
      pathFiltDir <- paste(deseq2FilesDirParent(), 
                        "/", 
                        input$deseq2FilesDirFiltPrefix, 
                        "_padj", 
                        padj_d(), 
                        "_log2FC", 
                        lfc_d(),
                        "/", sep="")
      
      # Create directories for new files
      dir.create(uniqueGeneListDir, recursive = TRUE)
      dir.create(pathFiltDir, recursive = TRUE)
      
      if (input$deseq2FilesRownames) {
        deseq2FilesRownames <- 1
      } else {
        deseq2FilesRownames <- NULL
      }
      
      filesSaved <- 0
      
      deseq2FileList <- list.files(pathRaw, full.names = TRUE)
      finalFileList <- deseq2FileList[which(basename(deseq2FileList) %in% input$deseq2FilesSelect)]
      
      for (rawFile in finalFileList) {
        
        fileExt <- tools::file_ext(rawFile)
        fileName <- tools::file_path_sans_ext(basename(rawFile))
        
        if (fileExt == "tsv") {
          fileSep <- "\t"
        } else if (fileExt == "csv") {
          fileSep <- ","
        }
        
        rawFile_df <- read.table(file = rawFile, 
          header=input$deseq2FilesHeader, 
          row.names = deseq2FilesRownames,
          sep=fileSep, 
          fill=FALSE)
        
        genesToKeep <- which(!is.na(rawFile_df[, input$deseq2padjCol]) 
                             & rawFile_df[, input$deseq2padjCol]<=padj_d() 
                             & abs(as.numeric(rawFile_df[, input$deseq2LfcCol]))>=lfc_d())
        filtFile_df <- rawFile_df[genesToKeep,]
        
        fileNameFilt <- paste(pathFiltDir, 
                              fileName, 
                              "_padj", 
                              padj_d(), 
                              "_log2FC", 
                              lfc_d(),
                              ".", fileExt,
                              sep="")
        
        # Save each as new file
        if (!file.exists(fileNameFilt) || input$overwriteDeseq2FilesFilt) {
          write.table(data.frame("gene_id"=rownames(filtFile_df),filtFile_df), 
                      file=fileNameFilt, 
                      sep=fileSep, 
                      row.names = FALSE, 
                      col.names = TRUE,
                      quote=FALSE)
          
          filesSaved <- (filesSaved + 1)
          
          message(paste("Saved:", fileNameFilt, sep=" "))
          message(paste("File: ", fileName, "_padj", padj_d(), "_log2FC", lfc_d(), ".", fileExt,
                        " contains ", length(genesToKeep)," DEGs", sep=""))
          message("")
        }
      }
      
      message(paste("Saved", filesSaved, "new filtered DESeq2 results files", sep=" "))
      message("")
      
      # Save list of unique genes
      uniqueDEG_df <- uniqueGeneListAuto()
      
      fileNameDEGs <- paste(uniqueGeneListDir, 
                            input$uniqueGeneListFilePrefix,
                            "_padj", 
                            padj_d(), 
                            "_log2FC", 
                            lfc_d(),
                            ".csv",
                            sep="")
      
      if (!file.exists(fileNameDEGs) || input$overwriteUniqueGeneListFile) {
        write.table(uniqueDEG_df, 
                    file=fileNameDEGs, 
                    sep=",", 
                    row.names = FALSE, 
                    col.names = TRUE,
                    quote=FALSE)
        message(paste("Saved:", fileNameDEGs, sep=" "))
        message(paste("File: ", input$uniqueGeneListFilePrefix, "_padj", padj_d(), "_log2FC", lfc_d(), ".csv",
                      " contains ", nrow(uniqueDEG_df)," unique DEGs", sep=""))
        message("")
      }
      
      update_saveFilteredData(update_saveFilteredData() + 1)
    })
    
    ########################################
    #     Heatmap - Server                 #
    ########################################
    
    uniqueTreatmentsWSamples_df <- reactive({
      req(input$metaData)
      
      # Import metadata
      metaData <- metaData()
      
      # Use treatment group mean gene count
      if (input$hm_useGroupMean){
        # Generate list of unique treatments to take average of
        uniqueTreatments <- expand.grid(lapply(metaData, unique))
        rownames(uniqueTreatments) <- apply(expand.grid(lapply(metaData, unique)), 1, paste, collapse = "_")
      
        # Identify which sample belongs to which treatment
        uniqueTreatments$samples <- ""
        for (i in 1:nrow(metaData)) {
          treatment <- paste(metaData[i,], collapse = "_")
          sample <- rownames(metaData)[i]
        
          if (uniqueTreatments[treatment, "samples"] == "") {
            uniqueTreatments[treatment, "samples"] <- sample
          } else {
            uniqueTreatments[treatment, "samples"] <- paste(uniqueTreatments[treatment, "samples"], sample, sep = ",")
          }
        }
      } else {
        # User per-sample gene count
        uniqueTreatments <- metaData
        rownames(uniqueTreatments) <- apply(cbind(sub(".*_", "", rownames(uniqueTreatments)), uniqueTreatments), 1, paste, collapse = "_")
        uniqueTreatments$samples <- rownames(metaData)
      }
      
      return(uniqueTreatments)
    })
    
    # Establish base requirements
    baseDataReady <- reactive ({
      req(input$countFile)
      req(input$metaData)
      req(isTruthy(input$deseq2FilesDir) || isTruthy(input$uniqueGeneListInput))
      
      TRUE
    })
    
    # Create raw count df for heatmap
    heatmap_df <- reactive({
      req(baseDataReady())
      
      geneCounts <- countData()
      
      
      if (isTruthy(input$deseq2FilesDir) && isTruthy(input$uniqueGeneListInput)) {
        if (input$hm_useGeneListAuto) {
          geneList <- uniqueGeneListAuto()
          message("Displaying heatmap using unique DEG list generated from input DESeq2 files, according to specified parameters.")
        } else {
          geneList <- uniqueGeneListInputTable()
          message(paste("Displaying heatmap using 'unique gene list' input file ", input$uniqueGeneListInput$name, ".", sep=""))
        }
      } else if (isTruthy(input$deseq2FilesDir)) {
        geneList <- uniqueGeneListAuto()
        message("Displaying heatmap using unique DEG list generated from input DESeq2 files, according to specified parameters. Previously generated gene list not found.")
      } else {
        geneList <- uniqueGeneListInputTable()
        message(paste("Displaying heatmap using 'unique gene list' input file ", input$uniqueGeneListInput$name, ". Input DESeq2 files not found.", sep=""))
      }
      
      
      uniqueTreatmentsWSamples <- uniqueTreatmentsWSamples_df()
      
      # Begin an empty dataframe with the list of desired genes
      combinedGeneCounts <- geneList
      rownames(combinedGeneCounts) = geneList$gene_id
      
      # Add a new column to the empty data frame for each sample name
      for (i in rownames(uniqueTreatmentsWSamples)) {
        combinedGeneCounts[,i] <- NA
        
        # Identify which samples belong to the new column
        samples <- str_split_1(uniqueTreatmentsWSamples[i, "samples"], ",")
        for (j in 1:length(samples)) {
          samples[j] <- make.names(samples[j])
        }
        
        # For each gene, take the mean gene count for the specified samples
        for (gene in rownames(combinedGeneCounts)) {
          combinedGeneCounts[gene, i] <- mean(as.numeric(geneCounts[gene, which(colnames(geneCounts) %in% samples)]))
        }
      }
      
      # Remove redundant gene_id column
      combinedGeneCounts$gene_id <- NULL
      
      return(combinedGeneCounts)
      
      #hm_data<-data()
      
      # Remove white space, replace empty gene names with NA
      #hm_data[, input$geneCol] <- gsub("^$", NA, trimws(hm_data[, input$geneCol]))
      
      # Find duplicate gene names
      #dup_geneNames <- unique(hm_data[duplicated(hm_data[, input$geneCol]), input$geneCol])
      
      # Replace NA and duplicate gene names with ID or name_ID, respectively
      #hm_data$names2use <- ifelse(is.na(hm_data[, input$geneCol]), hm_data[, input$geneAltCol], ifelse(hm_data[, input$geneCol] %in% dup_geneNames, paste0(hm_data[, input$geneCol], "_", hm_data[, input$geneAltCol]), hm_data[, input$geneCol]))
      
      # Replace rownames with unique, non-NA labels
      #rownames(hm_data) <- hm_data$names2use
      #hm_data$names2use <- NULL
      
      
      #if (input$lfc_direction == "both") {
      #  heatmap_tmp <- hm_data[which(!is.na(hm_data[, input$padjCol]) & hm_data[, input$padjCol]<padj_d()_hm & abs(as.numeric(hm_data[, input$lfcCol]))>=lfc_d()_hm),] 
      #} else if (input$lfc_direction == "negative") {
      #  heatmap_tmp <- hm_data[which(!is.na(hm_data[, input$padjCol]) & hm_data[, input$padjCol]<padj_d()_hm & as.numeric(hm_data[, input$lfcCol])<=-(lfc_d()_hm)),]
      #} else if (input$lfc_direction == "positive") {
      #  heatmap_tmp <- hm_data[which(!is.na(hm_data[, input$padjCol]) & hm_data[, input$padjCol]<padj_d()_hm & as.numeric(hm_data[, input$lfcCol])>=lfc_d()_hm),]
      #}
      
      #if (input$hm_sortBy == "padj") {
      #  heatmap_df <- heatmap_tmp[order(heatmap_tmp[,input$padjCol]),]
      #} else if (input$hm_sortBy == "log2FoldChange") {
      #  heatmap_df <- heatmap_tmp[order(heatmap_tmp[,input$lfcCol]),]
      #}
      
    })
    
    # Transform base df according to desired scaling
    scaled_heatmap_df <- reactive({
      req(heatmap_df())
      
      # Select desired genes to display from heatmap_df
      heatmap_df <- heatmap_df()
      #heatmap_df_selection <- heatmap_df_selection[1:input$hm_topXgenes,input$countCols]
      
      if (input$hm_scaling == "absolute TPM"){
        scaled_heatmap_df <- as.matrix(abs(heatmap_df))
      } else if (input$hm_scaling == "log2(TPM+1)") {
        # Calculate log2(mean+1) transformation of counts
        scaled_heatmap_df <- as.matrix(log2(heatmap_df+1))
      } else if (input$hm_scaling == "z-score") {
        # Calculate z-score transformation of counts
        t_heatmap_df <- t(as.matrix(heatmap_df))
        zt_heatmap_df <- scale(t_heatmap_df, center=T, scale = T)
        scaled_heatmap_df <- t(zt_heatmap_df)
      } else if (input$hm_scaling == "z-score(log2(TPM+1))") {
        log2_t_heatmap_df <- t(as.matrix(log2(heatmap_df+1)))
        log2_zt_heatmap_df <- scale(log2_t_heatmap_df, center=T, scale = T)
        scaled_heatmap_df <- t(log2_zt_heatmap_df)
      }
      
      
      return(scaled_heatmap_df)
    })
    
    # output$hm_topXgenes <- renderUI({
    #   req(heatmap_df())
    #   sliderInput("hm_topXgenes",
    #               "Choose number of genes to display:",
    #               min=1,
    #               max=nrow(heatmap_df()),
    #               value=nrow(heatmap_df()),
    #               step=1)
    # })
    
    # Divide columns based on clustering
    output$hm_cutreeCols <- renderUI({
      req(baseDataReady())
      
      choices=sequence(length(colnames(heatmap_df())), from=2L)
      
      selectizeInput("hm_cutreeCols",
                  "Split columns into X clusters:",
                  c("None", choices),
                  selected="None",
                  multiple=FALSE)
    })

    # Generate title based on input parameters
    output$hm_title_ui <- renderUI({
      req(baseDataReady())
      
      if (isTruthy(input$hm_showGeneClusters)){
        clusterTitle <- paste(" + Row clusters (k=", input$hm_geneClusterK, ")", sep="")
      } else {
        clusterTitle <- ""
      }
      
      n <- nrow(scaled_heatmap_df())
      
      if (input$hm_useGroupMean) {
        sampleValues <- "Group mean"
      } else {
        sampleValues <- "Per-sample"
      }
      
      textInput("hm_title","Diagram title:",
                value=paste("DEG heatmap (", sampleValues, ", ", input$hm_scaling, ")", clusterTitle, " (n=", n, ")", sep = ""))
    })
    
    # Create list of esthetic parameters for plotting
    hm_input_params <- reactive({
      list(
        color_rev    = input$hm_color_rev,
        color_number = input$hm_color_number,
        color        = input$hm_color,
        cluster_dist = input$hm_cluster_dist,
        gene_names   = input$hm_gene_names,
        sample_names = input$hm_sample_names,
        base_fontsize = input$hm_base_fontsize,
        gene_fontsize = input$hm_gene_fontsize,
        sample_fontsize = input$hm_sample_fontsize,
        sampleName_angle = input$hm_sampleName_angle
      )
    })
    
    # Create function to generate heatmap
    generate_pheatmap <- function(heatmap_df, metaData, dendrograms, title, colSplits, rowSplits, params, geneClusters = NA) {
      
      if(dendrograms == "both") {clusterCol=TRUE;clusterRow=TRUE}
      if(dendrograms == "none") {clusterCol=FALSE;clusterRow=FALSE}
      if(dendrograms == "column") {clusterCol=TRUE;clusterRow=FALSE}
      if(dendrograms == "row") {clusterCol=FALSE;clusterRow=TRUE}
      
      max_val <- max(heatmap_df, na.rm=TRUE)
      n_breaks <- 100
      brks <- (seq(0, 1, length.out = n_breaks)^2) * max_val
      
      if (params$color_rev) {
        usePalette <- colorRampPalette(rev(brewer.pal(n = params$color_number, name = params$color)))(100)
      } else {
        usePalette <- colorRampPalette(brewer.pal(n = params$color_number, name = params$color))(100)
      }
      
      if (rowSplits != "None"){
        if (clusterRow){
          cutTreeRows <- as.numeric(rowSplits)
        } else {
          cutTreeRows <- NA
          message("Please choose to cluster genes in order to show cluster separations")
        }
      } else {
        cutTreeRows <- NA
      }
      
      if (colSplits != "None"){
        if (clusterCol){
          cutTreeCols <- as.numeric(colSplits)
        } else {
          cutTreeCols <- NA
          message("Please choose to cluster samples in order to show cluster separations")
        }
      } else {
        cutTreeCols <- NA
      }
        

      
      p <- pheatmap(
        mat = heatmap_df,
        color = usePalette,
        breaks = brks,
        border_color = "grey60",
        cellwidth = NA,
        cellheight = NA,
        scale = "none",
        cluster_cols = clusterCol, 
        cluster_rows = clusterRow,
        clustering_distance_rows = params$cluster_dist,
        clustering_distance_cols = params$cluster_dist,
        clustering_method = "complete",
        cutree_rows = cutTreeRows,
        cutree_cols = cutTreeCols,
        legend = TRUE,
        annotation_row = geneClusters,
        annotation_col = metaData,
        annotation_legend = TRUE,
        annotation_names_row = TRUE,
        annotation_names_col = TRUE,
        drop_levels = TRUE,
        show_rownames = params$gene_names,
        show_colnames = params$sample_names,
        main = title,
        fontsize = params$base_fontsize,
        fontsize_row = params$gene_fontsize,
        fontsize_col = params$sample_fontsize,
        angle_col = params$sampleName_angle,
        display_numbers = FALSE,
        silent = FALSE,
        na_col = "#DDDDDD"
      )
      
      return(p)
    }
    
    # Generate base heatmap plot
    hm_plot_base <- reactive({
      req(scaled_heatmap_df())
      req(hm_input_params())
      req(input$hm_title)
      req(input$hm_cutreeRows)
      req(input$hm_cutreeCols)
      req(input$hm_showDendrogram)
      
      # Remove 'samples' column from metadata
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      
      generate_pheatmap(heatmap_df = scaled_heatmap_df(), 
                        metaData = metaData, 
                        dendrograms = input$hm_showDendrogram, 
                        colSplits = input$hm_cutreeCols, 
                        rowSplits = input$hm_cutreeRows, 
                        title = input$hm_title,
                        params = hm_input_params()
                        )
    })
    
    # From base plot, extract cluster IDs 
    hm_geneCluster_df <- reactive({
      req(hm_plot_base())
      req(input$hm_geneClusterK)
      
      ph <- hm_plot_base()
      
      if (input$hm_showDendrogram %in% c("both", "row")){
        k <- input$hm_geneClusterK
        row_cluster <- cutree(ph$tree_row, k = k)
      
        hm_geneCluster_df <- data.frame(
          cluster = factor(row_cluster),
          row.names = names(row_cluster)
        )
      
        return(hm_geneCluster_df)
      } else {
        return(NULL)
      }
    })
    
    # Re-plot; use cluster IDs if specified by user
    hm_plot <- reactive({
      req(hm_plot_base())
      #req(scaled_heatmap_df())
      #req(input$hm_title)
      #req(input$hm_cutreeRows)
      #req(input$hm_cutreeCols)
      #req(input$hm_showDendrogram)
      
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      # Input scaled df
      scaled_heatmap_df <- scaled_heatmap_df()
      
      if (isTruthy(input$hm_showGeneClusters) && isTruthy(hm_geneCluster_df())) {
        geneClusters <- hm_geneCluster_df()
      } else {
        geneClusters <- NULL
      }
      
      generate_pheatmap(heatmap_df = scaled_heatmap_df, 
                        metaData = metaData, 
                        dendrograms = input$hm_showDendrogram, 
                        colSplits = input$hm_cutreeCols, 
                        rowSplits = input$hm_cutreeRows, 
                        title = input$hm_title,
                        params = hm_input_params(),
                        geneClusters = geneClusters
      )
    })
    
    # Output final plot
    output$HeatmapPlot <- renderPlot({
      req(hm_plot())
      #req(input$metaData)
      #req(isTruthy(input$deseq2FilesDir) || isTruthy(input$uniqueGeneListInput))
      
      hm_plot()
    }, height=800)
    
    # Generate cluster silhouette plot
    output$hm_clusterSilPlot <- renderPlot({
      req(baseDataReady())
      
      df <- scaled_heatmap_df()
      
      row_dist <- dist(df)
      sil_width <- numeric(10)
      
      for(k in 2:10){
        pam_fit <- pam(row_dist, k = k)
        sil_width[k] <- pam_fit$silinfo$avg.width
      }

      plot(1:10, sil_width, type="b", xlab="Number of clusters", ylab="Average silhouette width")
    }) %>% 
    bindEvent(input$hm_generateClusterSil)
    
    clusterFilesDirParent <- reactive({
      req(deseq2FilesDirParent())
      req(input$uniqueGeneListDir)
      
      uniqueGeneListDir <- paste(deseq2FilesDirParent(), 
                                 "/", input$uniqueGeneListDir, "/", 
                                 sep="")
      
      return(uniqueGeneListDir)
    })
    
    clusterFilesDir <- reactive({
      req(clusterFilesDirParent())
      
      if (input$hm_useGroupMean) {
        countType <- "group-mean"
      } else {
        countType <- "per-sample"
      }
      
      clusterFilesDir <- paste(clusterFilesDirParent(), "cluster_DEGlists_", 
            countType, "_", input$hm_scaling, "_", input$hm_cluster_dist,
            "_padj", padj_d(), "_lfc", lfc_d(), "/", sep="")
      
      return(clusterFilesDir) 
    })
    
    output$geneClusterFilesPath <- renderPrint({
      if (isTruthy(input$deseq2FilesDirParent)){
        clusterFilesDir()
      } else {
        "Please choose a parent directory in the 'Data processing' tab"
      }
    })
    
    # If cluster file folder exists and contains files, display warning
    clusterFileList <- reactive({
      req(clusterFilesDir())
      update_saveMainClusters()
      
      path <- clusterFilesDir()
      
      list.files(path, full.names = FALSE)
    }) #%>% 
    #bindEvent(input$clusterFilesDirParent)
    
    output$clusterWarningOverwriteFiles <- renderText({
      req(clusterFilesDir())
      
      existingFileList <- clusterFileList()
      
      if (length(existingFileList) == 0) {
        return("Cluster file directory is currently empty")
      } else if (input$overwriteMainClusterFiles) {
        return("WARNING: Cluster file directory currently contains the following files. Identical filenames will be overwritten")
      } else {
        return("WARNING: Cluster file directory currently contains the following files. Identical filenames will be ignored")
      }
    })
    
    # Display current cluster files
    output$clusterMainExistingFiles <- renderPrint({
      req(clusterFilesDir())
      
      existingFileList <- clusterFileList()
      
      if (length(existingFileList) == 0) {
        return("No files to display")
      } else {
        return(cat(paste(existingFileList, collapse = "\n")))
      }
    })
    
    # Reactive value to force reactive updates on file save
    update_saveMainClusters <- reactiveVal(0)
    
    saveGeneClusterFiles <- observeEvent(input$hm_saveGeneClusterFiles, {
      req(isTruthy(hm_geneCluster_df()))
      req(clusterFilesDir())
      message("Saving main heatmap cluster-membership and per-cluster DEG files")
      message("")
      
      k <- input$hm_geneClusterK
      
      cluster_df <- hm_geneCluster_df()
      cluster_df <- data.frame("gene_id"=rownames(cluster_df),cluster_df)
      clusterFilesDir <- clusterFilesDir()
      
      dir.create(clusterFilesDir,showWarnings=FALSE, recursive=TRUE)
      
      fileNameCluster_df <- paste(clusterFilesDirParent(), 
                            input$uniqueGeneListFilePrefix,
                            "_cluster_membership_padj", 
                            padj_d(), 
                            "_log2FC", 
                            lfc_d(), 
                            "_k", 
                            k,
                            ".csv",
                            sep="")
      # Save full cluster gene list
      if (!file.exists(fileNameCluster_df) || input$overwriteMainClusterFiles) {
        write.table(cluster_df, 
                    file=fileNameCluster_df, 
                    sep=",", 
                    row.names = FALSE, 
                    col.names = TRUE,
                    quote=FALSE)
        message(paste("Saved main heatmap cluster-membership file:", fileNameCluster_df, sep=" "))
        message("")
      
      }
      # Save per-cluster gene lists
      for (i in 1:k) {
        genes_i <- cluster_df[which(cluster_df$cluster == i),]
        
        filename_i <- paste(clusterFilesDir,
                            input$hm_clusterFilePrefix,
                            "_m-k",
                            i, "of", k,
                            sep="")
        if (!file.exists(paste(filename_i, ".txt", sep="")) || input$overwriteMainClusterFiles) {
          if (input$saveMainClusterTxt) {
            writeLines(genes_i$gene_id, paste(filename_i, ".txt", sep=""))
          
            message(paste("Saved per-cluster DEG file: ", filename_i, ".txt", sep=""))
            message("")
          }
        }
        if (!file.exists(paste(filename_i, ".csv", sep="")) || input$overwriteMainClusterFiles) {
          write.table(genes_i,
                  file=paste(filename_i, ".csv", sep=""),
                  sep=",", 
                  row.names = FALSE, 
                  col.names = TRUE,
                  quote=FALSE)
          
          message(paste("Saved per-cluster DEG file: ", filename_i, ".csv", sep=""))
          message("")
        }
      }
      
      update_saveMainClusters(update_saveMainClusters() + 1)
    })
    
    output$hm_filename <- renderUI({
      req(baseDataReady())
      
      if (isTruthy(input$hm_showGeneClusters)){
        clusterTitle <- paste("-clustered_k", input$hm_geneClusterK, sep="")
      } else {
        clusterTitle <- ""
      }
      
      if (input$hm_useGroupMean) {
        sampleValues <- "group_mean"
      } else {
        sampleValues <- "per_sample"
      }
      
      textInput("hm_filename", "Filename", 
                value = paste("DEG_main_heatmap-", sampleValues, clusterTitle, sep = ""))
    })
    
    # Download heatmap as pdf
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0(input$hm_filename, ".", input$hm_file_format)
      },
      content = function(file) {
        width  <- input$hm_plot_width
        height <- input$hm_plot_height
        dpi    <- input$hm_plot_dpi
        
        switch(input$hm_file_format,
               "pdf" = {
                 pdf(file, width = width, height = height)
                 print(hm_plot())
                 dev.off()
               },
               "png" = {
                 png(file, width = width, height = height, units = "in", res = dpi)
                 print(hm_plot())
                 dev.off()
               },
               "svg" = {
                 svg(file, width = width, height = height)
                 print(hm_plot())
                 dev.off()
               },
               "tiff" = {
                 tiff(file, width = width, height = height, units = "in", res = dpi)
                 print(hm_plot())
                 dev.off()
               }
        )
      }
    )
    
    ########################################
    #     Sub-cluster heatmap - Server     #
    ########################################
    output$hmSUB_clusterChoice <- renderUI({
      selectizeInput("hmSUB_clusterChoice",
                     "Select main cluster # to display",
                     choices=c(sequence(input$hm_geneClusterK)),
                     multiple=FALSE)
    })
    
    # Order heatmap data according to main plot output
    ordered_heatmap_df <- reactive({
      req(hm_plot_base())
      
      scaled_heatmap_df <- scaled_heatmap_df()
      ph <- hm_plot()
      
      if (input$hm_showDendrogram %in% c("both", "row")){
        rowOrder <- rownames(scaled_heatmap_df)[ph$tree_row[["order"]]]
      } else {
        rowOrder <- rownames(scaled_heatmap_df)
      }
      if (input$hm_showDendrogram %in% c("both", "column")) { 
        colOrder <- colnames(scaled_heatmap_df)[ph$tree_col[["order"]]]
      } else {
        colOrder <- colnames(scaled_heatmap_df)
      }
      
      return(scaled_heatmap_df[rowOrder, colOrder])
    })
    
    # Keep only specified clusters
    reduced_heatmap_df <- reactive({
      req(ordered_heatmap_df())
      req(input$hmSUB_clusterChoice)
      req(isTruthy(hm_geneCluster_df()))
      
      ordered_heatmap_df <- ordered_heatmap_df()
      
      cluster_df <- hm_geneCluster_df()
      
      genesToKeep <- rownames(cluster_df)[which(cluster_df$cluster == as.numeric(input$hmSUB_clusterChoice))]
      
      return(ordered_heatmap_df[which(rownames(ordered_heatmap_df) %in% genesToKeep),])
    })
    
    output$hmSUB_cutreeCols <- renderUI({
      req(reduced_heatmap_df())
      
      choices=sequence(length(colnames(reduced_heatmap_df())), from=2L)
      
      selectizeInput("hmSUB_cutreeCols",
                     "Split columns into X clusters:",
                     c("None", choices),
                     selected="None",
                     multiple=FALSE)
    })
    
    output$hmSUB_title_ui <- renderUI({
      req(reduced_heatmap_df())
      req(input$hmSUB_clusterChoice)
      
      if (isTruthy(input$hmSUB_showGeneClusters)){
        clusterTitle <- paste(" + Row clusters (k=", input$hmSUB_geneClusterK, ")", sep="")
      } else {
        clusterTitle <- ""
      }
      
      n <- nrow(reduced_heatmap_df())
      
      if (input$hm_useGroupMean) {
        sampleValues <- "Group mean"
      } else {
        sampleValues <- "Per-sample"
      }
      
      textInput("hmSUB_title","Diagram title:",
                value=paste("DEG heatmap cluster ", input$hmSUB_clusterChoice, "/", input$hm_geneClusterK, " (", sampleValues, ", ", input$hm_scaling, ")", clusterTitle, " (n=", n, ")", sep = ""))
    })
    
    # Plot cluster-specific heatmap
    hmSUB_plot_base <- reactive({
      req(reduced_heatmap_df())
      req(hm_input_params())
      req(input$hmSUB_clusterChoice)
      req(input$hmSUB_title)
      req(input$hmSUB_cutreeRows)
      req(input$hmSUB_cutreeCols)
      
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      
      reduced_heatmap_df <- reduced_heatmap_df()
      
      
      generate_pheatmap(heatmap_df = reduced_heatmap_df, 
                        metaData = metaData, 
                        dendrograms = input$hmSUB_showDendrogram, 
                        colSplits = input$hmSUB_cutreeCols, 
                        rowSplits = input$hmSUB_cutreeRows, 
                        title = input$hmSUB_title,
                        params = hm_input_params()
                        )
    })
    
    hmSUB_geneCluster_df <- reactive({
      req(hmSUB_plot_base())
      req(input$hmSUB_geneClusterK)
      
      ph <- hmSUB_plot_base()
      
      if (input$hmSUB_showDendrogram %in% c("both", "row")){
        k <- input$hmSUB_geneClusterK
        row_cluster <- cutree(ph$tree_row, k = k)
        
        hmSUB_geneCluster_df <- data.frame(
          cluster = factor(row_cluster),
          row.names = names(row_cluster)
        )
        
        return(hmSUB_geneCluster_df)
      } else {
        return(NULL)
      }
      
    })
    
    hmSUB_plot <- reactive({
      req(hmSUB_plot_base())
      
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      
      # Input scaled df
      reduced_heatmap_df <- reduced_heatmap_df()
      
      if (isTruthy(input$hmSUB_showGeneClusters) && isTruthy(hmSUB_geneCluster_df())) {
        geneClusters <- hmSUB_geneCluster_df()
      } else {
        geneClusters <- NULL
      }
      
      generate_pheatmap(heatmap_df = reduced_heatmap_df,
                        metaData = metaData, 
                        dendrograms = input$hmSUB_showDendrogram, 
                        colSplits = input$hmSUB_cutreeCols, 
                        rowSplits = input$hmSUB_cutreeRows, 
                        title = input$hmSUB_title,
                        params = hm_input_params(),
                        geneClusters = geneClusters)
    })
    
    output$SUBHeatmapPlot <- renderPlot({
      
      hmSUB_plot()
    }, height=800)
  
    
    # Generate cluster silhouette plot
    output$hmSUB_clusterSilPlot <- renderPlot({
      req(reduced_heatmap_df())
      #req(input$countFile)
      #req(input$metaData)
      #req(isTruthy(input$deseq2FilesDir) || isTruthy(input$uniqueGeneListInput))
      
      df <- reduced_heatmap_df()
      
      row_dist <- dist(df)
      sil_width <- numeric(10)
      
      for(k in 2:10){
        pam_fit <- pam(row_dist, k = k)
        sil_width[k] <- pam_fit$silinfo$avg.width
      }
      
      plot(1:10, sil_width, type="b", xlab="Number of clusters", ylab="Average silhouette width")
    }) %>% 
    bindEvent(input$hmSUB_generateClusterSil)
    
    # Create SUBclusterFilesDir
    SUBclusterFilesDir <- reactive({
      req(clusterFilesDir())
      
      k <- input$hm_geneClusterK
      kc <- input$hmSUB_clusterChoice
      
      SUBclusterFilesDir <- paste(clusterFilesDir(), "sub-cluster_DEGlists_m-k", 
                                  kc, "of", k, "/", sep="")
      
      return(SUBclusterFilesDir)
    })
    
    output$geneSUBClusterFilesPath <- renderPrint({
      if (isTruthy(input$deseq2FilesDirParent)){
        SUBclusterFilesDir()
      } else {
        "Please choose a parent directory in the 'Data processing' tab"
      }
    })
    
    # If cluster file folder exists and contains files, display warning
    SUBclusterFileList <- reactive({
      req(SUBclusterFilesDir())
      update_saveSubClusters()
      
      path <- SUBclusterFilesDir()
      
      list.files(path, full.names = FALSE)
    }) #%>% 
    #bindEvent(input$clusterFilesDirParent)
    
    
    output$SUBclusterWarningOverwriteFiles <- renderText({
      req(SUBclusterFilesDir())
      
      existingFileList <- SUBclusterFileList()
      
      if (length(existingFileList) == 0) {
        return("Sub-cluster file directory is currently empty")
      } else if (input$overwriteSUBClusterFiles) {
        return("WARNING: Sub-cluster file directory currently contains the following files. Identical filenames will be overwritten")
      } else {
        return("WARNING: Sub-cluster file directory currently contains the following files. Identical filenames will be ignored")
      }
    })
    
    # Display current cluster files
    output$SUBclusterMainExistingFiles <- renderPrint({
      req(SUBclusterFilesDir())
      
      existingFileList <- SUBclusterFileList()
      
      if (length(existingFileList) == 0) {
        return("No files to display")
      } else {
        return(cat(paste(existingFileList, collapse = "\n")))
      }
    })
    
    # Reactive value to force reactive updates on file save
    update_saveSubClusters <- reactiveVal(0)
    
    saveGeneSUBClusterFiles <- observeEvent(input$hmSUB_saveGeneClusterFiles, {
      req(isTruthy(hmSUB_geneCluster_df()))
      req(SUBclusterFilesDir())
      message("Saving sub heatmap cluster-membership and per-cluster DEG files")
      message("")
      
      km <- input$hm_geneClusterK
      kc <- input$hmSUB_clusterChoice
      ksub <- input$hmSUB_geneClusterK
      
      subcluster_df <- hmSUB_geneCluster_df()
      subcluster_df <- data.frame("gene_id"=rownames(subcluster_df),subcluster_df)
      subclusterFilesDir <- SUBclusterFilesDir()
      
      dir.create(subclusterFilesDir,showWarnings=FALSE, recursive=TRUE)
      
      fileNamesubCluster_df <- paste(clusterFilesDir(), 
                                  input$hm_clusterFilePrefix,
                                  "_sub-cluster_membership_m-k",
                                  kc, "of", km,
                                  "_s-k", ksub,
                                  ".csv",
                                  sep="")
      
      # Save full cluster gene list
      if (!file.exists(fileNamesubCluster_df) || input$overwriteSUBClusterFiles) {
        write.table(subcluster_df, 
                    file=fileNamesubCluster_df, 
                    sep=",", 
                    row.names = FALSE, 
                    col.names = TRUE,
                    quote=FALSE)
        message(paste("Saved sub heatmap cluster-membership file:", fileNamesubCluster_df, sep=" "))
        message("")
        
      }
      # Save per-cluster gene lists
      for (i in 1:ksub) {
        genes_i <- subcluster_df[which(subcluster_df$cluster == i),]
        
        filename_i <- paste(subclusterFilesDir,
                            input$hmSUB_clusterFilePrefix,
                            "_m-k", 
                            kc, "of", km,
                            "_s-k",
                            i, "of", ksub,
                            sep="")
        if (!file.exists(paste(filename_i, ".txt", sep="")) || input$overwriteSUBClusterFiles) {
          if (input$saveSUBClusterTxt) {  
            writeLines(genes_i$gene_id, paste(filename_i, ".txt", sep=""))
          
            message(paste("Saved per-sub-cluster DEG file: ", filename_i, ".txt", sep=""))
            message("")
          }
        }
        if (!file.exists(paste(filename_i, ".csv", sep="")) || input$overwriteSUBClusterFiles) {
          write.table(genes_i,
                      file=paste(filename_i, ".csv", sep=""),
                      sep=",", 
                      row.names = FALSE, 
                      col.names = TRUE,
                      quote=FALSE)
          
          message(paste("Saved per-sub-cluster DEG file: ", filename_i, ".csv", sep=""))
          message("")
        }
      }
      
      update_saveSubClusters(update_saveSubClusters() + 1)
    })
    
    output$hmSUB_filename <- renderUI({
      req(reduced_heatmap_df())
      req(input$hmSUB_clusterChoice)
      
      if (isTruthy(input$hmSUB_showGeneClusters)){
        clusterTitle <- paste("-clustered_k", input$hmSUB_geneClusterK, sep="")
      } else {
        clusterTitle <- ""
      }
      
      if (input$hm_useGroupMean) {
        sampleValues <- "group_mean"
      } else {
        sampleValues <- "per_sample"
      }
      
      textInput("hmSUB_filename", "Filename", 
                value = paste("DEG_sub_heatmap-", sampleValues, "-main_cluster", input$hmSUB_clusterChoice, clusterTitle, sep = ""))
    })
    
    # Download heatmap as pdf
    output$download_SUBheatmap <- downloadHandler(
      filename = function() {
        paste0(input$hmSUB_filename, ".", input$hmSUB_file_format)
      },
      content = function(file) {
        width  <- input$hmSUB_plot_width
        height <- input$hmSUB_plot_height
        dpi    <- input$hmSUB_plot_dpi
        
        switch(input$hmSUB_file_format,
               "pdf" = {
                 pdf(file, width = width, height = height)
                 print(hmSUB_plot())
                 dev.off()
               },
               "png" = {
                 png(file, width = width, height = height, units = "in", res = dpi)
                 print(hmSUB_plot())
                 dev.off()
               },
               "svg" = {
                 svg(file, width = width, height = height)
                 print(hmSUB_plot())
                 dev.off()
               },
               "tiff" = {
                 tiff(file, width = width, height = height, units = "in", res = dpi)
                 print(hmSUB_plot())
                 dev.off()
               }
        )
      }
    )
    
    
    ########################################
    #     Sub-sub-cluster heatmap - Server #
    ########################################
    output$hmSUB2_clusterChoice <- renderUI({
      selectizeInput("hmSUB2_clusterChoice",
                     "Select sub cluster # to display",
                     choices=c(sequence(input$hmSUB_geneClusterK)),
                     multiple=FALSE)
    })
    
    # Order heatmap data according to main plot output
    ordered_SUBheatmap_df <- reactive({
      req(hmSUB_plot_base())
      
      reduced_heatmap_df <- reduced_heatmap_df()
      ph <- hmSUB_plot()
      
      if (input$hmSUB_showDendrogram %in% c("both", "row")){
        rowOrder <- rownames(reduced_heatmap_df)[ph$tree_row[["order"]]]
      } else {
        rowOrder <- rownames(reduced_heatmap_df)
      }
      if (input$hmSUB_showDendrogram %in% c("both", "column")) { 
        colOrder <- colnames(reduced_heatmap_df)[ph$tree_col[["order"]]]
      } else {
        colOrder <- colnames(reduced_heatmap_df)
      }
      
      return(reduced_heatmap_df[rowOrder, colOrder])
    })
    
    # Keep only specified clusters
    reduced_SUBheatmap_df <- reactive({
      req(ordered_SUBheatmap_df())
      req(input$hmSUB2_clusterChoice)
      req(isTruthy(hmSUB_geneCluster_df()))
      
      ordered_SUBheatmap_df <- ordered_SUBheatmap_df()
      
      cluster_df <- hmSUB_geneCluster_df()
      
      genesToKeep <- rownames(cluster_df)[which(cluster_df$cluster == as.numeric(input$hmSUB2_clusterChoice))]
      
      return(ordered_SUBheatmap_df[which(rownames(ordered_SUBheatmap_df) %in% genesToKeep),])
    })
    
    output$hmSUB2_cutreeCols <- renderUI({
      req(reduced_SUBheatmap_df())
      #req(input$countFile)
      #req(input$metaData)
      #req(isTruthy(input$deseq2FilesDir) || isTruthy(input$uniqueGeneListInput))
      
      choices=sequence(length(colnames(reduced_SUBheatmap_df())), from=2L)
      
      selectizeInput("hmSUB2_cutreeCols",
                     "Split columns into X clusters:",
                     c("None", choices),
                     selected="None",
                     multiple=FALSE)
    })
    
    output$hmSUB2_title_ui <- renderUI({
      req(reduced_SUBheatmap_df())
      req(input$hmSUB2_clusterChoice)
      
      #req(**Use group mean or individual samples)
      
      if (isTruthy(input$hmSUB2_showGeneClusters)){
        clusterTitle <- paste(" + Row clusters (k=", input$hmSUB2_geneClusterK, ")", sep="")
      } else {
        clusterTitle <- ""
      }
      
      n <- nrow(reduced_SUBheatmap_df())
      
      if (input$hm_useGroupMean) {
        sampleValues <- "Group mean"
      } else {
        sampleValues <- "Per-sample"
      }
      
      textInput("hmSUB2_title","Diagram title:",
                value=paste("DEG heatmap main cluster ", input$hmSUB_clusterChoice, "/", input$hm_geneClusterK,  
                            "; sub-cluster ", input$hmSUB2_clusterChoice, "/", input$hmSUB_geneClusterK,
                            " (", sampleValues, ", ", input$hm_scaling, ")",
                            clusterTitle, " (n=", n, ")", sep = ""))
    })
    
    # Plot cluster-specific heatmap
    hmSUB2_plot_base <- reactive({
      req(reduced_SUBheatmap_df())
      req(hm_input_params())
      req(input$hmSUB2_clusterChoice)
      req(input$hmSUB2_title)
      req(input$hmSUB2_cutreeRows)
      req(input$hmSUB2_cutreeCols)
      
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      
      reduced_SUBheatmap_df <- reduced_SUBheatmap_df()
      
      
      generate_pheatmap(heatmap_df = reduced_SUBheatmap_df, 
                        metaData = metaData, 
                        dendrograms = input$hmSUB2_showDendrogram, 
                        colSplits = input$hmSUB2_cutreeCols, 
                        rowSplits = input$hmSUB2_cutreeRows, 
                        title = input$hmSUB2_title,
                        params = hm_input_params()
      )
    })
    
    hmSUB2_geneCluster_df <- reactive({
      req(hmSUB2_plot_base())
      req(input$hmSUB2_geneClusterK)
      
      ph <- hmSUB2_plot_base()
      
      if (input$hmSUB2_showDendrogram %in% c("both", "row")){
        k <- input$hmSUB2_geneClusterK
        row_cluster <- cutree(ph$tree_row, k = k)
        
        hmSUB2_geneCluster_df <- data.frame(
          cluster = factor(row_cluster),
          row.names = names(row_cluster)
        )
        
        return(hmSUB2_geneCluster_df)
      } else {
        return(NULL)
      }
      
    })
    
    hmSUB2_plot <- reactive({
      req(hmSUB2_plot_base())
      #req(reduced_heatmap_df())
      #req(input$hmSUB_clusterChoice)
      #req(input$hmSUB_title)
      #req(input$hmSUB_cutreeRows)
      #req(input$hmSUB_cutreeCols)
      
      metaData <- uniqueTreatmentsWSamples_df()
      metaData$samples <- NULL
      # Input scaled df
      reduced_SUBheatmap_df <- reduced_SUBheatmap_df()
      
      if (isTruthy(input$hmSUB2_showGeneClusters) && isTruthy(hmSUB2_geneCluster_df())) {
        geneClusters <- hmSUB2_geneCluster_df()
      } else {
        geneClusters <- NULL
      }
      
      generate_pheatmap(heatmap_df = reduced_SUBheatmap_df,
                        metaData = metaData, 
                        dendrograms = input$hmSUB2_showDendrogram, 
                        colSplits = input$hmSUB2_cutreeCols, 
                        rowSplits = input$hmSUB2_cutreeRows, 
                        title = input$hmSUB2_title,
                        params = hm_input_params(),
                        geneClusters = geneClusters)
    })
    
    output$SUB2HeatmapPlot <- renderPlot({
      
      hmSUB2_plot()
    }, height=800)
    
    # Generate cluster silhouette plot
    output$hmSUB2_clusterSilPlot <- renderPlot({
      req(reduced_SUBheatmap_df())
      
      df <- reduced_SUBheatmap_df()
      
      row_dist <- dist(df)
      sil_width <- numeric(10)
      
      for(k in 2:10){
        pam_fit <- pam(row_dist, k = k)
        sil_width[k] <- pam_fit$silinfo$avg.width
      }
      
      plot(1:10, sil_width, type="b", xlab="Number of clusters", ylab="Average silhouette width")
    }) %>% 
      bindEvent(input$hmSUB2_generateClusterSil)
    
    # Create SUB2clusterFilesDir
    SUB2clusterFilesDir <- reactive({
      req(SUBclusterFilesDir())
      
      k <- input$hm_geneClusterK
      kc <- input$hmSUB_clusterChoice
      k2 <- input$hmSUB_geneClusterK
      k2c <- input$hmSUB_clusterChoice
      
      SUB2clusterFilesDir <- paste(SUBclusterFilesDir(), "sub2-cluster_DEGlists_m-k", 
                                  kc, "of", k, 
                                  "_sub-k", k2c, "of", k2,"/", sep="")
      
      return(SUB2clusterFilesDir)
    })
    
    output$geneSUB2ClusterFilesPath <- renderPrint({
      if (isTruthy(input$deseq2FilesDirParent)){
        SUB2clusterFilesDir()
      } else {
        "Please choose a parent directory in the 'Data processing' tab"
      }
    })
    
    # If cluster file folder exists and contains files, display warning
    SUB2clusterFileList <- reactive({
      req(SUB2clusterFilesDir())
      update_saveSub2Clusters()
      
      path <- SUB2clusterFilesDir()
      
      list.files(path, full.names = FALSE)
    }) #%>% 
      #bindEvent(input$clusterFilesDirParent)
    
    
    output$SUB2clusterWarningOverwriteFiles <- renderText({
      req(SUB2clusterFilesDir())
      
      existingFileList <- SUB2clusterFileList()
      
      if (length(existingFileList) == 0) {
        return("Sub2-cluster file directory is currently empty")
      } else if (input$overwriteSUB2ClusterFiles) {
        return("WARNING: Sub2-cluster file directory currently contains the following files. Identical filenames will be overwritten")
      } else {
        return("WARNING: Sub2-cluster file directory currently contains the following files. Identical filenames will be ignored")
      }
    })
    
    # Display current cluster files
    output$SUB2clusterMainExistingFiles <- renderPrint({
      req(SUB2clusterFilesDir())
      
      existingFileList <- SUB2clusterFileList()
      
      if (length(existingFileList) == 0) {
        return("No files to display")
      } else {
        return(cat(paste(existingFileList, collapse = "\n")))
      }
    })
    
    # Reactive value to force reactive updates on file save
    update_saveSub2Clusters <- reactiveVal(0)
    
    saveGeneSUB2ClusterFiles <- observeEvent(input$hmSUB2_saveGeneClusterFiles, {
      req(isTruthy(hmSUB2_geneCluster_df()))
      req(SUB2clusterFilesDir())
      message("Saving sub2 heatmap cluster-membership and per-cluster DEG files")
      message("")
      
      km <- input$hm_geneClusterK
      kc <- input$hmSUB_clusterChoice
      
      ksub <- input$hmSUB_geneClusterK
      ksubc <- input$hmSUB2_clusterChoice
      
      ksub2 <- input$hmSUB2_geneClusterK
      
      subcluster_df <- hmSUB2_geneCluster_df()
      subcluster_df <- data.frame("gene_id"=rownames(subcluster_df),subcluster_df)
      subclusterFilesDir <- SUB2clusterFilesDir()
      
      dir.create(subclusterFilesDir,showWarnings=FALSE, recursive=TRUE)
      
      fileNamesubCluster_df <- paste(SUBclusterFilesDir(), 
                                     input$hmSUB_clusterFilePrefix,
                                     "_sub-cluster_membership_m-k", 
                                     kc, "of", km,
                                     "_s-k", ksubc, "of", ksub,
                                     "_s2-k", ksub2,
                                     ".csv",
                                     sep="")
      
      # Save full cluster gene list
      if (!file.exists(fileNamesubCluster_df) || input$overwriteSUBClusterFiles) {
        write.table(subcluster_df, 
                    file=fileNamesubCluster_df, 
                    sep=",", 
                    row.names = FALSE, 
                    col.names = TRUE,
                    quote=FALSE)
        message(paste("Saved sub2 heatmap cluster-membership file:", fileNamesubCluster_df, sep=" "))
        message("")
        
      }
      
      # Save per-cluster gene lists
      for (i in 1:ksub2) {
        genes_i <- subcluster_df[which(subcluster_df$cluster == i),]
        
        filename_i <- paste(subclusterFilesDir,
                            input$hmSUB2_clusterFilePrefix,
                            "_m-k", kc, "of", km,
                            "_s-k_", ksubc, "of", ksub,
                            "_sub2-k", i, "of", ksub2,
                            sep="")
        if (!file.exists(paste(filename_i, ".txt", sep="")) || input$overwriteSUBClusterFiles) {
          if (input$saveSUB2ClusterTxt){
            writeLines(genes_i$gene_id, paste(filename_i, ".txt", sep=""))
            
            message(paste("Saved per-sub2-cluster DEG file: ", filename_i, ".txt", sep=""))
            message("")
          }
        }
        if (!file.exists(paste(filename_i, ".csv", sep="")) || input$overwriteSUBClusterFiles) {
          write.table(genes_i,
                      file=paste(filename_i, ".csv", sep=""),
                      sep=",", 
                      row.names = FALSE, 
                      col.names = TRUE,
                      quote=FALSE)
          
          message(paste("Saved per-sub2-cluster DEG file: ", filename_i, ".csv", sep=""))
          message("")
        }
      }
      
      update_saveSub2Clusters(update_saveSub2Clusters() + 1)
    })
    
    output$hmSUB2_filename <- renderUI({
      req(reduced_SUBheatmap_df())
      req(input$hmSUB2_clusterChoice)
      
      if (isTruthy(input$hmSUB2_showGeneClusters)){
        clusterTitle <- paste("-clustered_k", input$hmSUB2_geneClusterK, sep="")
      } else {
        clusterTitle <- ""
      }
      
      if (input$hm_useGroupMean) {
        sampleValues <- "group_mean"
      } else {
        sampleValues <- "per_sample"
      }
      
      textInput("hmSUB2_filename", "Filename", 
                value = paste("DEG_sub2_heatmap-", sampleValues, 
                              "-main_k", input$hmSUB_clusterChoice, 
                              "-sub_k", input$hmSUB2_clusterChoice, 
                              clusterTitle, sep = ""))
    })
    
    # Download heatmap as pdf
    output$download_SUB2heatmap <- downloadHandler(
      filename = function() {
        paste0(input$hmSUB2_filename, ".", input$hmSUB2_file_format)
      },
      content = function(file) {
        width  <- input$hmSUB2_plot_width
        height <- input$hmSUB2_plot_height
        dpi    <- input$hmSUB2_plot_dpi
        
        switch(input$hmSUB2_file_format,
               "pdf" = {
                 pdf(file, width = width, height = height)
                 print(hmSUB2_plot())
                 dev.off()
               },
               "png" = {
                 png(file, width = width, height = height, units = "in", res = dpi)
                 print(hmSUB2_plot())
                 dev.off()
               },
               "svg" = {
                 svg(file, width = width, height = height)
                 print(hmSUB2_plot())
                 dev.off()
               },
               "tiff" = {
                 tiff(file, width = width, height = height, units = "in", res = dpi)
                 print(hmSUB2_plot())
                 dev.off()
               }
        )
      }
    )
    
    
    ########################################
    #     Cluster GO-terms - Server        #
    ########################################
    # Input GO <-> gene conversion df
    # Validate file ending
    go2geneExtension <- reactive({
      req(input$go2geneConversion)
      go2gene <- input$go2geneConversion
      go2geneExt <- tools::file_ext(go2gene$datapath)
      
      validate(need(go2geneExt %in% c("tsv", "csv", "txt"), "Please upload either a csv, tsv, or txt file"))
      
      return(go2geneExt)
    })
    
    output$go2geneTxtExt <- renderUI({
      req(go2geneExtension())

      go2geneExt <- go2geneExtension()
      
      if (go2geneExt == "txt") {
        selectizeInput("go2geneTxtExt",
                       "Select text file delimiter:",
                       choices=c("tab"="\t",
                                 "comma"=","),
                       multiple=FALSE)
      }
      
    })
    
    go2geneDelim <- reactive({
      req(go2geneExtension() %in% c("tsv","csv") || isTruthy(input$go2geneTxtExt))
      
      if (go2geneExtension() %in% c("tsv","csv")){
        if (go2geneExtension() == "tsv") {
          go2geneDelim <- "\t"
        } else if (go2geneExtension() == "csv") {
          go2geneDelim <- ","
        }
      } else {
        go2geneDelim <- input$go2geneTxtExt
      }
      
      return(go2geneDelim)
    })
    
    go2gene <- reactive({
      req(go2geneDelim())
      
      go2gene <- input$go2geneConversion
      
      go2geneSep <- go2geneDelim()
      
      if (input$go2geneRownames) {
        go2geneRownames <- 1
      } else {
        go2geneRownames <- NULL
      }
      
      read.csv(
        go2gene$datapath, 
        header=input$go2geneHeader, 
        row.names = go2geneRownames,
        sep=go2geneSep, 
        fill=FALSE)
    })  
    
    # Input GO ID <-> term conversion df
    # Validate file ending
    go2termExtension <- reactive({
      req(input$go_termConversion)
      go2term <- input$go_termConversion
      go2termExt <- tools::file_ext(go2term$datapath)
      
      validate(need(go2termExt %in% c("tsv", "csv", "txt"), "Please upload either a csv, tsv, or txt file"))
      
      return(go2termExt)
    })
    
    output$go2termTxtExt <- renderUI({
      req(go2termExtension())
      
      go2termExt <- go2termExtension()
      
      if (go2termExt == "txt") {
        selectizeInput("go2termTxtExt",
                       "Select text file delimiter:",
                       choices=c("tab"="\t",
                                 "comma"=","),
                       multiple=FALSE)
      }
      
    })
    
    go2termDelim <- reactive({
      req(go2termExtension() %in% c("tsv","csv") || isTruthy(input$go2termTxtExt))
      
      if (go2termExtension() %in% c("tsv","csv")){
        if (go2termExtension() == "tsv") {
          go2termDelim <- "\t"
        } else if (go2termExtension() == "csv") {
          go2termDelim <- ","
        }
      } else {
        go2termDelim <- input$go2termTxtExt
      }
      
      return(go2termDelim)
    })
    
    go2term <- reactive({
      req(go2termDelim())
      
      go2term <- input$go_termConversion
      
      go2termSep <- go2termDelim()
      
      if (input$go2termRownames) {
        go2termRownames <- 1
      } else {
        go2termRownames <- NULL
      }
      
      read.csv(
        go2term$datapath, 
        header=input$go2termHeader, 
        row.names = go2termRownames,
        sep=go2termSep)
    }) 
    
    
    # Select appropriate total gene list based on heatmap choice
    go_baseGeneList <- reactive ({
      if (input$go_hmChoice == "Main"){
        go_base_df <- scaled_heatmap_df()
      } else if (input$go_hmChoice == "Sub"){
        go_base_df <- reduced_heatmap_df()
      } else if (input$go_hmChoice == "Sub2"){
        go_base_df <- reduced_SUBheatmap_df()
      }
      
      return(go_base_df)
    })
    
    # Select appropriate cluster df based on heatmap choice
    go_baseClusterDF <- reactive ({
      if (input$go_hmChoice == "Main"){
        go_cluster_df <- hm_geneCluster_df()
      } else if (input$go_hmChoice == "Sub"){
        go_cluster_df <- hmSUB_geneCluster_df()
      } else if (input$go_hmChoice == "Sub2"){
        go_cluster_df <- hmSUB2_geneCluster_df()
      }
      
      return(go_cluster_df)
    })
    
    # Identify desired cluster(s) from heatmap choice
    output$go_clusterChoice <- renderUI({
      if (input$go_hmChoice == "Main"){
        clust_choices <- c(sequence(input$hm_geneClusterK))
      } else if (input$go_hmChoice == "Sub"){
        clust_choices <- c(sequence(input$hmSUB_geneClusterK))
      } else if (input$go_hmChoice == "Sub2"){
        clust_choices <- c(sequence(input$hmSUB2_geneClusterK))
      }
      
      selectizeInput("go_clusterChoice",
                     "Choose cluster # from relevant heatmap",
                     choices=clust_choices,
                     selected=clust_choices[1],
                     multiple=TRUE)
    })
    
    # Process choices to create gene list according to chosen heatmap and cluster(s)
    go_refinedGeneList <- reactive({
      req(input$go_clusterChoice)
      
      gene_df <- go_baseGeneList()
      
      cluster_df <- go_baseClusterDF()
      
      genesToKeep <- rownames(cluster_df)[which(cluster_df$cluster %in% input$go_clusterChoice)]
      
      return(gene_df[which(rownames(gene_df) %in% genesToKeep),])
    })
    
    # Refine GO:IDs based on desired genes
    go_refinedIDList <- reactive({
      req(input$go2geneConversion)
      req(go_refinedGeneList())
      
      go2gene <- go2gene()
      
      go_refinedIDList <- go2gene[which(go2gene[,2] %in% rownames(go_refinedGeneList())),]
      
      return(go_refinedIDList)
    })
    
    go_EGO <- reactive({
      req(go_refinedIDList())
      req(input$go_termConversion)
      
      go_refinedList <- go_refinedIDList()
      go2gene <- go2gene()
      go2term <- go2term()
      
      cluster_geneList <- unique(go_refinedList[,2])
      all_genes <- unique(go2gene[,2]) 
      
      ego <- enricher(
        gene = cluster_geneList,                          ## Which genes found in both specified gene cluster AND go2gene list
        universe = all_genes,                             ## Unique gene names in go2gene
        TERM2GENE = go2gene,                              ## go2gene
        TERM2NAME = go2term,                              ## go2term
        pAdjustMethod = input$ego_pAdjustMethod,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = input$ego_minGSSize,
        maxGSSize = input$ego_maxGSSize
      )
      
      return(ego)
    })
    
    go_EGO_df <- reactive({
      req(go_EGO())
      
      ego_df <- as.data.frame(go_EGO())
      ego_df$GeneRatio_num <- sapply(strsplit(ego_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      ego_df$neglogFDR <- -log10(ego_df$p.adjust)
      
      return(ego_df)
    })
    
    # Decide number of GO terms to include in 'preview' plot
    output$go_previewLimit <- renderUI({
      req(go_EGO_df())
      
      go_EGO_df <- go_EGO_df()
      
      value <- min(20, nrow(go_EGO_df))
      
      numericInput("go_previewLimit",
                   "# terms to display; Full plot",
                   min = 1,
                   max = nrow(go_EGO_df),
                   value = value,
                   step = 1)
    })
    
    go_previewImage <- reactive({
      req(input$go_previewLimit)
      req(input$go_largePlotStyle)
      
      ego <- go_EGO()
      
      topN <- input$go_previewLimit
      
      if (input$go_largePlotStyle == "bubble"){
        return(dotplot(ego, showCategory = topN) + ggtitle(input$go_largePlotTitle))
      } else if (input$go_largePlotStyle == "bar"){
        return(barplot(ego, showCategory = topN) + ggtitle(input$go_largePlotTitle))
      }
    })
    
    output$go_previewPlot <- renderPlot({
      
      go_previewImage()
    })
    
    output$goPreview_filename <- renderUI({
      req(input$go_previewLimit)
      req(input$go_largePlotStyle)
      
      textInput("goPreview_filename", "Filename", 
                value = paste("GO_term-preview_", input$go_largePlotStyle, 
                              "_plot-top_", input$go_previewLimit, sep = ""))
    })
    
    # Download heatmap as pdf
    output$download_goPreview <- downloadHandler(
      filename = function() {
        paste0(input$goPreview_filename, ".", input$goPreview_file_format)
      },
      content = function(file) {
        width  <- input$goPreview_plot_width
        height <- input$goPreview_plot_height
        dpi    <- input$goPreview_plot_dpi
        
        switch(input$goPreview_file_format,
               "pdf" = {
                 pdf(file, width = width, height = height)
                 print(go_previewImage())
                 dev.off()
               },
               "png" = {
                 png(file, width = width, height = height, units = "in", res = dpi)
                 print(go_previewImage())
                 dev.off()
               },
               "svg" = {
                 svg(file, width = width, height = height)
                 print(go_previewImage())
                 dev.off()
               },
               "tiff" = {
                 tiff(file, width = width, height = height, units = "in", res = dpi)
                 print(go_previewImage())
                 dev.off()
               }
        )
      }
    )
    
    # Decide count limit to display on histogram
    output$go_countLimit <- renderUI({
      req(go_EGO_df())
      
      go_EGO_df <- go_EGO_df()
      
      value <- min(10, max(go_EGO_df$Count))
      
      numericInput("go_countLimit",
                   "Minimum count number to plot",
                   min = 1,
                   max = max(go_EGO_df$Count),
                   value = value,
                   step = 1)
    })
    
    # Display histogram of gene counts
    output$go_countHistogram <- renderPlot({
      req(input$go_countLimit)
      
      go_EGO_df <- go_EGO_df()
      
      xmax <- round_any(max(go_EGO_df$Count), 10, ceiling)
      ticks <- xmax/10
      
      hist((go_EGO_df$Count)[which(go_EGO_df$Count > input$go_countLimit)],
           breaks=(max(go_EGO_df$Count)/2),
           main="Frequency of GO terms with a given gene count",
           xlab="Gene count",
           ylab="# of GO terms",
           xlim=c(0,xmax),
           xaxp=c(0,xmax,ticks)
      )
    })
    
    ########################################
    #     Cluster bubble plot - Server     #
    ########################################
    # Small bubble plot df preparation
    go_EGO_plotDF <- reactive({
      req(go_EGO_df())
      
      ego_df <- go_EGO_df()
      
      # Filter by desired p and q values
      ego_df <- ego_df[which(ego_df$p.adjust < input$ego_pvalueCutoff & ego_df$qvalue < input$ego_qvalueCutoff),]
      
      # Order df according to choice
      ego_df_ordered <- ego_df[order(ego_df[,input$ego_dfSort], decreasing = input$ego_dfSortDesc),]
      
      return(ego_df_ordered)
    })
    
    # Decide number of GO terms to include in bubble plot
    output$go_bubbleLimit <- renderUI({
      req(go_EGO_plotDF())
      
      go_EGO_df <- go_EGO_plotDF()
      
      value <- min(15, nrow(go_EGO_df))
      
      numericInput("go_bubbleLimit",
                   "Top X terms to display (from arranged GO table)",
                   min = 1,
                   max = nrow(go_EGO_df),
                   value = value,
                   step = 1)
    })
    
    # Bubble plot manual selection
    output$go_bubbleTermChoice <- renderUI({
      req(go_EGO_plotDF())
      
      if (input$go_bubbleChooseTerms){
        df <- go_EGO_plotDF()
        
        choices = df$Description
      
        selectizeInput("go_bubbleTermChoice", "Choose GO terms to display:", choices=choices, selected = NULL, multiple=TRUE)
      }
    })
    
    # Bubble plot title
    output$go_bubbleTitle <- renderUI({
      req(input$go_bubbleLimit)
      
      if (input$go_bubbleChooseTerms) {
        textInput("go_bubbleTitle","Diagram title:",
                  value="GO enrichment (Selected terms)")
      } else {
        textInput("go_bubbleTitle","Diagram title:",
                  value=paste("GO enrichment (Top ", input$go_bubbleLimit, ")", sep = ""))
      }
    })
    
    # Reorder bubble plot df
    go_EGO_plotDFordered <- reactive({
      req(input$input$go_bubbleLimit)
      req(!input$go_bubbleChooseTerms || length(input$go_bubbleTermChoice) > 0)
            
      ego_df <- go_EGO_plotDF()
      
      # Filter to # of terms
      # if (input$go_bubbleChooseTerms) {
      #   ego_plot <- ego_df[which(ego_df$Description %in% input$go_bubbleChooseTerms),]
      # } else {
      #   
      # }
      
      ego_plot <- ego_df[1:input$go_bubbleLimit,]
      
      ego_plot <- ego_plot[order(ego_plot[,input$ego_plotSort], decreasing = input$ego_plotSortDesc),]
      ego_plot$ID <- factor(ego_plot$ID, levels = rev(ego_plot$ID))  
      
      return(ego_plot)
    })
    
    # Small bubble plot
    go_EGO_plotImage <- reactive({
      req(go_EGO_plotDFordered())
      req(input$go_bubbleTitle)
      
      ego_plot <- go_EGO_plotDFordered()
      label_map <- setNames(ego_plot$Description, ego_plot$ID)
      
      midpoint <- max(ego_plot[,input$ego_dfSort])/2
      
      colourLegend <- switch(input$ego_dfSort,
                             RichFactor = "Rich Factor",
                             FoldEnrichment = "Fold Enrichment",
                             GeneRatio_num = "Gene ratio",
                             neglogFDR = expression(-log[10] ~ p.adjust),
                             input$ego_dfSort  # default
      )
      
      p <- ggplot(ego_plot, aes(x = GeneRatio_num, y = ID))
      
      
      p <- p + switch(input$ego_bubblePlotTheme,
                      bw = theme_bw(base_size = input$ego_themeBaseSize),
                      linedraw = theme_linedraw(base_size = input$ego_themeBaseSize),
                      light = theme_light(base_size = input$ego_themeBaseSize),
                      dark = theme_dark(base_size = input$ego_themeBaseSize),
                      minimal = theme_minimal(base_size = input$ego_themeBaseSize),
                      classic = theme_classic(base_size = input$ego_themeBaseSize),
                      gray = theme_gray(base_size = input$ego_themeBaseSize)  # default / "gray"
      )
      
      p <- p + 
        theme(
          legend.position = input$ego_bubbleLegendPosition,
          legend.box.spacing = unit(input$ego_legendSpacing, "pt"),
          legend.margin = margin_auto(input$ego_legendMargin),
          legend.box = input$ego_bubbleLegendStacking,
          legend.title = element_text(size = input$ego_legendTitleSize),
          legend.text = element_text(size = input$ego_legendTextSize),
          plot.title = element_text(size = input$ego_titleSize),
          axis.text.y = element_text(size = input$ego_yAxisTextSize),
          axis.text.x = element_text(size = input$ego_xAxisTextSize),
          plot.margin = margin(t = 5, r = 30, b = 5, l = 5)
        )
      
      p <- p +
        scale_x_continuous(expand = expansion(mult = c(0.09, 0.03))) +
        scale_y_discrete(labels = ~ str_wrap(label_map[.x], input$ego_yAxisStrWrap)) +
        scale_size_continuous(name = "Gene count") +
        labs(title = input$go_bubbleTitle) +
        coord_cartesian(clip = "off")
      
      if (input$ego_revPlotColours) {
        colour_low <- input$go_bubbleColorHigh
        colour_high <- input$go_bubbleColorLow
      } else {
        colour_low <- input$go_bubbleColorLow
        colour_high <- input$go_bubbleColorHigh
      }
      
      if (input$go_bubbleDotOutline) {
        p <- p + 
          geom_point(aes(size = Count, fill = .data[[input$ego_dfSort]]), 
                     alpha = input$ego_geomAlpha, 
                     shape = 21, 
                     stroke = input$ego_geomOutlineSize, 
                     colour = "black") +
          scale_fill_gradient2(
            name = colourLegend,
            low  = colour_low,
            mid = input$go_bubbleColorMid,
            high = colour_high,
            midpoint = midpoint,
            limits = c(0, NA)
          )
      } else {
        p <- p + 
          geom_point(aes(size = Count, colour = .data[[input$ego_dfSort]]), 
                     alpha = input$ego_geomAlpha) +
          scale_colour_gradient2(
            name = colourLegend,
            low  = colour_low,
            mid = input$go_bubbleColorMid,
            high = colour_high,
            midpoint = midpoint,
            limits = c(0, NA)
          )
      }
      
      p <- p + labs(
        y = if (input$go_yAxisShowTitle) input$go_yAxisTitle else NULL,
        x = if (input$go_xAxisShowTitle) input$go_xAxisTitle else NULL
      )
      
      p <- p + theme(
        axis.ticks.y       = if (input$go_yAxisTicks)    element_line()  else element_blank(),
        axis.ticks.x       = if (input$go_xAxisTicks)    element_line()  else element_blank(),
        panel.grid.major.y = if (input$go_yAxisGridlines) element_line() else element_blank(),
        panel.grid.major.x = if (input$go_xAxisGridlines) element_line() else element_blank()
      )

      if (input$ego_legendBox) {
        p <- p + theme(legend.box.background = element_rect(), 
                       legend.background = element_blank(),
                       legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0))
      }
      
      inside_theme <- switch(input$ego_bubbleLegendPositionInside,
                             "bottom left" = theme(
                               legend.justification.inside = c("left",  "bottom"),
                               legend.position.inside = c(0, 0)
                             ),
                             "bottom right" = theme(
                               legend.justification.inside = c("right", "bottom"),
                               legend.position.inside = c(1, 0)
                             ),
                             "top left" = theme(
                               legend.justification.inside = c("left",  "top"),
                               legend.position.inside = c(0, 1)
                             ),
                             "top right" = theme(
                               legend.justification.inside = c("right", "top"),
                               legend.position.inside = c(1, 1)
                             )
      )
      
      if (input$ego_bubbleLegendPosition == "inside") {
        p <- p + inside_theme
      } else {
        p <- p + theme(legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0))
      }
      
      return(p)
    })
    
    output$go_EGO_plot <- renderPlot({
      
      go_EGO_plotImage()
    })
    
    output$goEGO_filename <- renderUI({
      req(input$go_bubbleLimit)
      
      textInput("goEGO_filename", "Filename", 
                value = paste("GO_term-bubble_plot-top_", input$go_bubbleLimit,
                              "-", input$ego_plotSort, sep = ""))
    })
    
    # Download bubble plot
    output$download_goEGO <- downloadHandler(
      filename = function() {
        paste0(input$goEGO_filename, ".", input$goEGO_file_format)
      },
      content = function(file) {
        width  <- input$goEGO_plot_width
        height <- input$goEGO_plot_height
        dpi    <- input$goEGO_plot_dpi
        
        switch(input$goEGO_file_format,
               "pdf" = {
                 pdf(file, width = width, height = height)
                 print(go_EGO_plotImage())
                 dev.off()
               },
               "png" = {
                 png(file, width = width, height = height, units = "in", res = dpi)
                 print(go_EGO_plotImage())
                 dev.off()
               },
               "svg" = {
                 svg(file, width = width, height = height)
                 print(go_EGO_plotImage())
                 dev.off()
               },
               "tiff" = {
                 tiff(file, width = width, height = height, units = "in", res = dpi)
                 print(go_EGO_plotImage())
                 dev.off()
               }
        )
      }
    )
  }
)
