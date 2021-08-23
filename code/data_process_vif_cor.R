library(spatstat)
library(here)
library(sp)
library(rgeos)
library(maptools)
library(GISTools)
library(tmap)
library(sf)
library(geojson)
library(geojsonio)
library(tmaptools)
library(tidyverse)
library(ggthemes)
library(mapview)
library(RColorBrewer)
library(spdep)
library(geojsonsf)
library(ggplot2)
library(ggalt)
library(ggsn)
library(spatialreg)
library(spgwr)
library(broom)
library(car)
library(scales)
library(BAMMtools)
library(reshape2)


library(corrplot)
library(mma)



df <- st_read(here::here("thesis","886_for_r_0809.csv"))

summary(df)

df$nearest_Length <- as.numeric(as.character(df$nearest_Length))
df$abc_mean <- as.numeric(as.character(df$abc_mean))
df$Point_Count <- as.numeric(as.character(df$Point_Count))
df$energy_tot <- as.numeric(as.character(df$energy_tot))
df$h_nutrients_calories_norm <- as.numeric(as.character(df$h_nutrients_calories_norm))
df$representativeness_norm <- as.numeric(as.character(df$representativeness_norm))
df$Obesity2006_2008 <- as.numeric(as.character(df$Obesity2006_2008))





#use the symbox() function in the car package to try a range of transfomations along Tukey’s ladder:
symbox(~nearest_Length, 
       df, 
       na.rm=T,
       powers=seq(-3,3,by=.5))

ggplot(df, aes(x=(nearest_Length)^0.5)) + 
  geom_histogram()
#power0.5 for proximity leads to normal distribution

symbox(~abc_mean, 
       df, 
       na.rm=T,
       powers=seq(-3,3,by=.5))

ggplot(df, aes(x=log(abc_mean))) + 
  geom_histogram()
#log for diversity leads to normal distribution

symbox(~Point_Count, 
       df, 
       na.rm=T,
       powers=seq(-3,3,by=.5))

ggplot(df, aes(x=(Point_Count)^0.5)) + 
  geom_histogram()
#transformation for diversity doesn't success


#apply transformation to data
df2 = select(df, 3:8)

summary(df2)

df2$nearest_Length_0_5 <- as.numeric((df$nearest_Length)^0.5)
df2$abc_mean_log <- as.numeric(log(df$abc_mean))



df3 = select(df2, 1:3, 6:8)

summary(df3)


#normalize

df3$energy_tot_norm <- rescale(df3$energy_tot, to = c(0, 1), from = range(df3$energy_tot, na.rm = TRUE, finite = TRUE)) 
df3$Obesity2006_2008_norm <- rescale(df3$Obesity2006_2008, to = c(0, 1), from = range(df3$Obesity2006_2008, na.rm = TRUE, finite = TRUE)) 
df3$Point_Count_norm <- rescale(df3$Point_Count, to = c(0, 1), from = range(df3$Point_Count, na.rm = TRUE, finite = TRUE)) 
df3$nearest_Length_0_5_norm <- rescale(df3$nearest_Length_0_5, to = c(0, 1), from = range(df3$nearest_Length_0_5, na.rm = TRUE, finite = TRUE)) 
df3$abc_mean_log_norm <- rescale(df3$abc_mean_log, to = c(0, 1), from = range(df3$abc_mean_log, na.rm = TRUE, finite = TRUE)) 



df4 = select(df3, 2, 7:11)
summary(df4)

#check histograms

gghist <- ggplot(df4, 
                 aes(x=energy_tot_norm)) + 
  geom_histogram(color="white", fill="#2b8cbe", alpha=0.5)+
  labs(title="Histogram of Energy", 
       x="Value", 
       y="Frequency")
# add a vertical line to the hisogram showing mean tempearture
gghist + geom_vline(aes(xintercept=mean(energy_tot_norm, 
                                        na.rm=TRUE)),
                    color="#2b8cbe", 
                    linetype="dashed", 
                    size=1)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


gghist <- ggplot(df4, 
                 aes(x=Point_Count_norm)) + 
  geom_histogram(color="white", fill="#2b8cbe", alpha=0.5)+
  labs(title="Histogram of Density", 
       x="Value", 
       y="Frequency")
# add a vertical line to the hisogram showing mean tempearture
gghist + geom_vline(aes(xintercept=mean(Point_Count_norm, 
                                        na.rm=TRUE)),
                    color="#2b8cbe", 
                    linetype="dashed", 
                    size=1)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

gghist <- ggplot(df4, 
                 aes(x=Obesity2006_2008_norm)) + 
  geom_histogram(color="white", fill="#2b8cbe", alpha=0.5)+
  labs(title="Histogram of Obesity", 
       x="Value", 
       y="Frequency")
# add a vertical line to the hisogram showing mean tempearture
gghist + geom_vline(aes(xintercept=mean(Obesity2006_2008_norm, 
                                        na.rm=TRUE)),
                    color="#2b8cbe", 
                    linetype="dashed", 
                    size=1)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))




#vif and correlation

#df4 is the final dataframe ready for analysis

#continuous y
#arrange column order 
df4$Obesity2006_2008_norm2 <- df4$Obesity2006_2008_norm
df5 = select(df4, 1:2, 4:7)

colnames(df5) <- c("Nutrient_entropy", "Energy", "Distance_power_0.5","Homogeneity_log","Density","Obesity_rate")

#correlation 
M2 = cor(df5)
corrplot(M2,
         tl.col = "black", method = 'color',addCoef.col = 'white',order = 'AOE')

#vif
model2 <- lm(Obesity2006_2008_norm2 ~ nearest_Length_0_5_norm + 
               abc_mean_log_norm +
               Point_Count_norm +
               h_nutrients_calories_norm +
               energy_tot_norm, data = df5)

vif(model2)




#re-selecting data to map variable distribution geographically
df7 = select(df, 1:8)
df7$nearest_Length_0_5 <- as.numeric((df$nearest_Length)^0.5)
df7$abc_mean_log <- as.numeric(log(df$abc_mean))
summary(df8)
df8 = select(df7, 1:5, 8:10)


df8$energy_tot_norm <- rescale(df8$energy_tot, to = c(0, 1), from = range(df8$energy_tot, na.rm = TRUE, finite = TRUE)) 
df8$Obesity2006_2008_norm <- rescale(df8$Obesity2006_2008, to = c(0, 1), from = range(df8$Obesity2006_2008, na.rm = TRUE, finite = TRUE)) 
df8$Point_Count_norm <- rescale(df8$Point_Count, to = c(0, 1), from = range(df8$Point_Count, na.rm = TRUE, finite = TRUE)) 
df8$nearest_Length_0_5_norm <- rescale(df8$nearest_Length_0_5, to = c(0, 1), from = range(df8$nearest_Length_0_5, na.rm = TRUE, finite = TRUE)) 
df8$abc_mean_log_norm <- rescale(df8$abc_mean_log, to = c(0, 1), from = range(df8$abc_mean_log, na.rm = TRUE, finite = TRUE)) 
summary(df8)
df9 = select(df8, 1:2, 4, 9:13)




london_pop <- st_read(here::here("thesis", "msoa_popcen_2011.gpkg"))
london <- st_read(here::here("thesis", "msoa_2011.gpkg"))

tm_shape(london)+ 
  tm_graticules(col = "grey", alpha=0.2)+
  tm_polygons("RGN11NM", palette = c( "#80cdc1", "#018571"), alpha = 0.45)+
  tm_layout(main.title = "London MSOA Map", legend.outside = TRUE)+
  tm_compass(north=0, size = 2)+
  tm_scale_bar(text.size=0.3)

london2 = select(london, 1:2, 13)
london_data <- left_join(london2, df9, by = c("MSOA11CD" = "area_id"))
london_data2 <- drop_na(london_data)


breaks2 = c(0, 0.2, 0.5, 0.6, 0.8, 1.0)



tm4 <- tm_shape(london_data2) +
  tm_polygons("energy_tot_norm",
              breaks = breaks2,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(b) Total Energy", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("energy_tot_norm",
              breaks = breaks2,
              palette="PuBu",
              title = "Energy Consumption in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t2=tmap_arrange(tm4, legend, ncol = 3)

t2


breaks3 = c(0.6, 0.625, 0.65, 0.675,0.7, 0.725)
tm3 <- tm_shape(london_data2) +
  tm_polygons("h_nutrients_calories_norm",
              breaks = breaks3,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(a) Nutrient Entropy", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("h_nutrients_calories_norm",
              breaks = breaks3,
              palette="PuBu",
              title = "Nutrient Entropy in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t3=tmap_arrange(tm3, legend, ncol = 3)
t3


breaks2 = c(0, 0.2, 0.5, 0.6, 0.8, 1.0)
tm5 <- tm_shape(london_data2) +
  tm_polygons("Obesity2006_2008_norm",
              breaks = breaks2,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(a) Obesity Rate", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("Obesity2006_2008_norm",
              breaks = breaks2,
              palette="PuBu",
              title = "Obesity in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t5=tmap_arrange(tm5, legend, ncol = 3)
t5


breaks2 = c(0, 0.2, 0.5, 0.6, 0.8, 1.0)
tm6 <- tm_shape(london_data2) +
  tm_polygons("Point_Count_norm",
              breaks = breaks2,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(a) Retail Density", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("Point_Count_norm",
              breaks = breaks2,
              palette="PuBu",
              title = "Retail Density in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t6=tmap_arrange(tm6, legend, ncol = 3)
t6

tm7 <- tm_shape(london_data2) +
  tm_polygons("nearest_Length_0_5_norm",
              breaks = breaks2,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(b) Retail Proximity", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("nearest_Length_0_5_norm",
              breaks = breaks2,
              palette="PuBu",
              title = "Retail Proximity in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t7=tmap_arrange(tm7, legend, ncol = 3)
t7

tm8 <- tm_shape(london_data2) +
  tm_polygons("abc_mean_log_norm",
              breaks = breaks2,
              palette="PuBu" )+
  tm_layout(frame=FALSE)+
  tm_layout(legend.show = FALSE)+
  tm_credits("(c) Retail Diversity", position=c(0,0.7), size=1.5)

legend <- tm_shape(london_data2) +
  tm_polygons("abc_mean_log_norm",
              breaks = breaks2,
              palette="PuBu",
              title = "Retail Diversity in London") +
  tm_scale_bar(position=c(0.2,0.04), text.size=0.6)+
  tm_compass(north=0, position=c(0.25,0.6))+
  tm_layout(legend.only = TRUE, legend.position=c(0.2,0.25),asp=0.1)


t8=tmap_arrange(tm8, legend, ncol = 3)
t8


#other analysis were done with SPSS and ArcGIS





