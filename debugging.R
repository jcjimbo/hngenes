

# Debugging Test for Mouse gene discrepancy 
# Barchart values don't align with the scatter plot for mouse 

# File used for barplot
b <- mouse_expression_data %>% filter(pair_id == "6295_ID")

# File used for scatter plot 
s <- m_expression_data2 %>% filter(pair_id == "6295_ID")

# scanning the app shows that the numbers between the two files are the same 
# appears to be an issue with the axis, specifically, with the barplot 
# No number is greater than 7! 




m_expression_data2 %>% filter(pair_id %in% "6295_ID") %>% 
    ggplot(aes(tissue, lognorm, fill = tissue)) + 
    geom_bar(stat = "summary") + ylab("log_normalised_counts") + xlab("") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    facet_grid(GeneID ~.)

# STAT SUMMARY!! 