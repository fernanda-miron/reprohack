## Cargamos las librerias que vamos a usar
library("vroom")
library("ggplot2")
library("dplyr")
library("tidyr")
library("gridExtra")

## Vamos a empezar con violin plots
## Leemos la data que vamos a utilizar
violin_data <- vroom("violin_plots.csv", col_names = T)

## Si observamos la figura del paper, solo estamos ploteando controles vs DLB
## Nuestro df tiene otra categoria llamada LBD, vamos a borrarla
## filter nos ayuda a filtrar, | significa OR
violin_data_clean <- violin_data %>% 
  filter(sampleType == "Control" | sampleType == "DLB")

## Una vez que tenemos los datos listos, vamos a graficar un violin plot 
## Vamos a empezar con cg04866173
## vamos a usar ggplot2
p1 <- ggplot(violin_data_clean, aes(x=as.factor(sampleType), y=cg04866173, fill=sampleType)) +
  geom_violin(width=1.4) + ## que tan ancha queremos la linea
  geom_boxplot(width=0.1, color="black", fill = "white") + ## que tan ancha queremos la linea
  scale_fill_manual(values = c("#6AB1DA","#EF6C68")) + ## que colores queremos poner ahi
  theme_bw() +
  theme(
    legend.position="none", ## borramos barra con legend
    plot.title = element_text(size=11) ## cambiamos tamanio de texto
  ) +
  ggtitle("A Violin wrapping a boxplot") + ## anadimos un titulo
  xlab("") ## borramos titulo en eje x
p1

## Pero la imagen final tiene multiples violin plots
## vamos a crear una funcion para automatizar los violin plots
plot_violin_cpg <- function(violin_data_clean, cpgs) { ## the function takes the dataframe and a vector
  
  plots <- list() ## abrimos una lista vacia
  
  for (cpg in cpgs) { ## dentro de la funcion creamos un loop for
    p <- ggplot(violin_data_clean, aes(x = as.factor(sampleType), y = .data[[cpg]], fill = sampleType)) +
      geom_violin(width = 1.4) +
      geom_boxplot(width = 0.1, color = "black", fill = "white") +
      scale_fill_manual(values = c("#6AB1DA", "#EF6C68")) +
      theme_light() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 11)
      ) +
      ggtitle(cpg) +
      xlab("") + 
      ylab("")
    
    plots[[cpg]] <- p ## guardamos el plot en la lista que empieza vacia e iteramos a traves de los cpgs
  }
  
  return(plots) ## lanzamos plots para poder tener la variable dentro del ambiente
}

## Vamos a usar la funcion
## Primero creamos un vector con los cpgs
cpgs <- c("cg16086807", "cg18800161", "cg16250093", "cg04866173", "cg11099930", "cg06951630", "cg24435966")

## Creamos los plots para cada cpg
plots <- plot_violin_cpg(violin_data_clean, cpgs)

## Usamos la funcion grid.arrange del paquete grid extra
## grob == graphical objects
final_plot <- grid.arrange(grobs = plots, ncol = 3)

## Finalmente guardamos
ggsave(filename = "boxplots.png", plot = final_plot, 
       height = 8, width = 12, bg = "white")

## Primero, vamos a leer los datos en R con vroom
methylation_data <- vroom("ontology.csv", col_names = F)

## Elegimos colores para la figura
color_map <- c("GO:BP" = "#FFAB0D", "GO:CC" = "#6AB1DA", "TF" = "#EF6C68")

## Making Figure
plot2 <- ggplot(methylation_data, aes(x = -log10(X3), y = X8, fill = as.factor(X7), size = X6)) +
  geom_point(shape = 21, color = "gray") +  # Change shape to 21 (a filled circle)
  scale_fill_manual(values = color_map) +
  scale_size_continuous(range = c(3, 10)) +
  scale_y_discrete(position = "right") +
  facet_wrap(~ X7, scales = "free_y", ncol = 1, strip.position = "left") +
  theme_light() +
  labs(
    x = expression(-log[10](italic(p))),
    y = NULL,
    size = "N. of genes",
    fill = "Source"  # Change color label to fill to match aesthetic
  ) +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_text(hjust = 1)
  )
plot2

## Finalmente guardamos
ggsave(filename = "ontology.png", plot = plot2, 
       height = 8, width = 10, bg = "white")


############# FINISH TEST
violin_data_group <- violin_data %>%
  select(cg04866173, sampleType) %>% 
  group_by(sampleType) %>% 
  summarise(mean = mean(cg04866173), probe = "cg04866173") %>% 
  pivot_wider(names_from = sampleType, values_from = mean) %>% 
  mutate(db_control_lbd = LBD - Control)