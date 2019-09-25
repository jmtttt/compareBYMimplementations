## code to prepare `DATASET` dataset goes here
shape_GER <- wrangle_shapefile(dir = "I:/Arbeitsgruppen/Survival/Daten/Shapefiles", layer = "VG250_KRS")
population <- wrangle_population("I:/Arbeitsgruppen/Survival/Daten/Populationsdaten/new/population_export.csv") %>% group_by(region_ID) %>% summarize(n = sum(bevoelkerung))

shape_GER@data <- dplyr::left_join(x = shape_GER@data,
                                  y = population,
                                  by = "region_ID")

#usethis::use_data("shape_GER")

save(shape_GER, file = "I:/Arbeitsgruppen/Survival/Nebenprojekte/MCMC vs. INLA/Simulation/Skripte/compareBYMimplementations/data/shape_GER.Rdata")
