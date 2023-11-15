library(tidyverse)
library(dplyr)
library(plotly)
library(countrycode)

# Download data
plasmodium <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Plasmodium&format=tsv")
anopheles <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Anopheles&format=tsv")
malaria <- read_csv("malaria_distribution.csv")

# Filter malaria causing Plasmodium species, and only obtain column of interest, Then obtain frequencies by country
plasmodium.Frequency <- plasmodium[plasmodium$species_name %in% c("Plasmodium falciparum", "Plasmodium vivax", "Plasmodium malariae", "Plasmodium ovale", "Plasmodium knowlesi"), c("species_name", "country")] %>%
    filter(!is.na(country)) %>%
    filter(!is.na(species_name)) %>%
    group_by(country) %>%
    count(country)

anopheles.Frequency <- anopheles[, c("species_name", "country")] %>%
    filter(!is.na(country)) %>%
    filter(!is.na(species_name)) %>%
    group_by(country) %>%
    count(country)

malaria.Frequency <- malaria[malaria$IsLatestYear == TRUE, c("Location", "FactValueNumeric")]

## Compiled this to immediately show cases of malaria with more than 0 as the fact value numeric
malaria.Frequency2 <- malaria[malaria$IsLatestYear == TRUE, c("Location", "FactValueNumeric")] %>% 
  + filter(FactValueNumeric >0)


# Change column names
colnames(plasmodium.Frequency) <- c("country", "plasmodium")
colnames(anopheles.Frequency) <- c("country", "anopheles")
colnames(malaria.Frequency) <- c("country", "malaria")

# Join all dataframes into a grand data frame where each country has malaria occurrence, sampling frequency of plasmodium and anopheles
total.Frequency <- left_join(x=malaria.Frequency, 
                             y=full_join(x=plasmodium.Frequency, 
                                         y=anopheles.Frequency, 
                                         by="country"), 
                             by="country") %>% filter(malaria > 0)

# Find continent where country is from
total.Frequency$continent <- countrycode(sourcevar=total.Frequency$country,
                                         origin="country.name",
                                         destination="continent")

# Dataframe of what percent of each continent is sampled, then create nested piechart
total.Frequency.Continent <- total.Frequency %>%
    group_by(continent) %>%
    summarise(malaria=sum(malaria, na.rm=TRUE),
              plasmodium=sum(plasmodium, na.rm=TRUE),
              anopheles=sum(anopheles, na.rm=TRUE),
              .groups='drop') %>%
    as.data.frame()

plot_ly(total.Frequency.Continent) %>%
    add_pie(labels=~continent, values=~malaria, type='pie', hole=0.660, sort=FALSE) %>%
    add_pie(total.Frequency.Continent, labels=~continent, values=~plasmodium, type='pie', hole=0.505, sort=FALSE,
            domain=list(x=c(0.175, 0.825), y=c(0.175, 0.825))) %>%
    add_pie(total.Frequency.Continent, labels=~continent, values=~anopheles, sort=FALSE,
            domain=list(x=c(0.34, 0.66), y=c(0.34, 0.66)))


##Create the data frame of frequencies from each genus and their respective continents. 
library(ggplot2)
Continent_data <- data.frame(
  continent = c("Africa","Americas","Asia", "Oceania"),
  malaria = c(163875665, 524158, 1247499, 736424),
  plasmodium = c(178, 190, 375, 278),
  anopheles = c(1813, 3847, 2916, 584)
)

# Reshaping the data for ggplot2 to display values of each category. I created a grouped bar chart to depict the data in an alternative way which may be more legible for some audiences. 
library(tidyr)
Continent_data_long <- gather(Continent_data, key = "variable", value = "value", -continent)

# Create the grouped bar chart. I used log to better visualize the data so it was more understandable and easier to read.
ggplot(Continent_data_long, aes(x = continent, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Grouped Bar Chart with Logarithmic Scale of total frequencies by each continent",
       x = "Continent", y = "Frequency Count", fill = "Genus of Mosquito") +
  scale_y_log10() +
  scale_fill_manual(values = c("malaria" = "blue", "plasmodium" = "black", "anopheles" = "purple")) +
  theme_minimal()

#Show unsampled of plasmodium and anopheles per continent and create pie chart. Edited so that there is a title along for the graph and a legend. Changed the trace titles to depict Anopheles vs Plasmodium. 
total.Count.Continent <- total.Frequency %>%
  group_by(continent) %>%
  summarise(unsampledPlasmodium=sum(is.na(plasmodium)),
            unsampledAnopheles=sum(is.na(anopheles)),
            .groups='drop') %>%
  as.data.frame() %>%
  filter(unsampledPlasmodium > 0) %>%
  filter(unsampledAnopheles > 0)

plot_ly(total.Count.Continent) %>%
  add_pie(labels=~continent, values=~unsampledPlasmodium, type='pie', hole=0.505, sort=FALSE, name = "Plasmodium") %>%
  add_pie(total.Count.Continent, labels=~continent, values=~unsampledAnopheles, 
          domain=list(x=c(0.25, 0.75), y=c(0.25, 0.75)),
          sort=FALSE, name = "Anopheles") %>%
  layout(
    title = "Nested Pie Chart of unsampled plasmodium and anopheles per continent",
    showlegend = TRUE,
    legend = list(
      x = 7,
      y = 0.5,
      traceorder = "normal",
      orientation = "v",
      title = list( 
        text = "Continents",
        side = "top")))


# Determine relationship between malaria occurrences and number of plasmodium and anopheles sampled
plasmodium.Frequency.Filtered <- total.Frequency[!is.na(total.Frequency$plasmodium), c("country", "continent", "malaria", "plasmodium")]
result <- cor.test(plasmodium.Frequency.Filtered$plasmodium,
                   plasmodium.Frequency.Filtered$malaria,
                   method="pearson")

result

anopheles.Frequency.Filtered <- total.Frequency[!is.na(total.Frequency$anopheles), c("country", "continent", "malaria", "anopheles")]
result <- cor.test(anopheles.Frequency.Filtered$anopheles,
                   anopheles.Frequency.Filtered$malaria,
                   method="pearson")

result

