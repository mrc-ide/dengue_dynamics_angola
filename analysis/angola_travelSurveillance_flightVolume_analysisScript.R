# Loading required libraries
library(dplyr); library(countrycode); library(ggplot2); library(tidyverse)

# Loading in pre-processed IATA data and filtering by Angola
iata_data <- readRDS("Q:/NAO_modelling/IATA_Data/yearly_IATA_data/iata_country_country_all_years.rds")
angola <- iata_data %>%
  filter(orig_country == "Angola" | dest_country == "Angola") %>%
  mutate(orig_simple = case_when(orig_country == "Portugal" ~ "Portugal",
                                 orig_country == "China" ~ "China",
                                 TRUE ~ "Other")) %>%
  mutate(dest_simple = case_when(dest_country == "Portugal" ~ "Portugal",
                                 dest_country == "China" ~ "China",
                                 TRUE ~ "Other")) 

# Calculating total volume of passengers arriving into Angola by year
angola_yearly_in <- angola %>%
  filter(dest_country == "Angola") %>%
  group_by(orig_simple, year) %>%
  summarise(passengers = sum(passengers))
angola_yearly_total_in_plot <- ggplot(angola_yearly_in, aes(x = year, y = passengers, col = orig_simple)) +
  geom_line() +
  lims(y = c(0, NA)) +
  labs(y = "Total Volume of Passengers Into Angola 2012-2021", x = "Year") +
  theme_bw()
ggsave(filename = "angola_yearly_total_in_plot.pdf", plot = angola_yearly_total_in_plot,
       width = 10, height = 5)

# Calculating total volume of passengers leaving Angola by year
angola_yearly_out <- angola %>%
  filter(orig_country == "Angola") %>%
  filter(dest_country != "Angola") %>%
  group_by(dest_simple, year) %>%
  summarise(passengers = sum(passengers))
angola_yearly_total_out_plot <- ggplot(angola_yearly_out, aes(x = year, y = passengers, col = dest_simple)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#E6491E", "#948D9B", "#459024")) +
  lims(y = c(0, NA)) +
  labs(y = "Total Volume of International\nDepartures from Angola 2012-2021", x = "") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(filename = "angola_yearly_total_out_plot.pdf", plot = angola_yearly_total_out_plot,
       width = 10, height = 5)

# Calculating the above, but by country
total_into_angola <- angola %>%
  filter(dest_country == "Angola") %>%
  group_by(orig_country) %>%
  summarise(passengers = sum(passengers)) %>%
  dplyr::arrange(desc(passengers))
x <- total_into_angola[-1, ]
x$passengers[1] / sum(x$passengers)
sum(total_into_angola[2:11, 2])

total_into_angola <- total_into_angola[2:11, ]
total_into_angola$iso <- countrycode::countrycode(total_into_angola$orig_country, origin = "country.name", destination = "iso3c")
total_into_angola$orig_country[5] <- "United Arab\nEmirates"
total_into_angola$orig_country[7] <- "United\nKingdom"
angola_yearly_total_in_plot_country <- ggplot(total_into_angola, aes(x = reorder(orig_country, passengers), y = passengers, fill = orig_country)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Total Volume of Passengers Into Angola 2012-2021", x = "") +
  theme(legend.position = "none") 
ggsave(filename = "angola_yearly_total_in_plot_country.pdf", plot = angola_yearly_total_in_plot_country,
       width = 7.5, height = 5)

total_out_angola <- angola %>%
  filter(orig_country == "Angola") %>%
  group_by(dest_country) %>%
  summarise(passengers = sum(passengers)) %>%
  dplyr::arrange(desc(passengers))
total_out_angola <- total_out_angola[2:11, ]
total_out_angola$iso <- countrycode::countrycode(total_out_angola$dest_country, origin = "country.name", destination = "iso3c")
total_out_angola$dest_country[4] <- "United Arab\nEmirates"
total_out_angola$dest_country[6] <- "United\nKingdom"
angola_yearly_total_out_plot_country <- ggplot(total_out_angola, aes(x = reorder(dest_country, passengers), y = passengers, fill = dest_country)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#948D9B",
                               "#E6491E",
                               "#948D9B",
                               "#948D9B",
                               "#948D9B",
                               "#948D9B",
                               "#459024",
                               "#948D9B",
                               "#948D9B",
                               "#948D9B"),
                    name = "Dest. Country") +
  labs(y = "Total Number of Arrivals\nfrom Angola 2012-2021", x = "",
       fill = "Dest. Country") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
sum(total_out_angola[2:11, 2])

flight_plot <- cowplot::plot_grid(angola_yearly_total_out_plot,
                                  angola_yearly_total_out_plot_country,
                                  labels = c("C", "D"), rel_widths = c(0.95, 1.25),
                                  align = "h", axis = "b")

ggsave(filename = "angola_yearly_total_out_plot_country.pdf", plot = angola_yearly_total_out_plot_country,
       width = 7.5, height = 5)

x <- total_out_angola[-1, ]
x$passengers[1] / sum(x$passengers)

cowplot::plot_grid(angola_yearly_total_out_plot_country, angola_yearly_total_in_plot_country,
                   labels = c("A", "B"))

## plotting geosentinel data
geosentinel <- read.csv("C:/Users/cw1716/Documents/Q_Drive_Copy2/Active_Research_Projects/angola_dengue/GeoSentinel ECCo Data Request_DENV Portugal_CEVDI-INSA (2).csv")

df <- geosentinel %>%
  filter(source == "original_Geosentinel") %>%
  select(Birth_Country, Country_of_Residence, Recent_Travel_Countries, Clinic_Visit_Date, Recent_Travel_Countries) %>%
  mutate(date = as.Date(Clinic_Visit_Date, format = "%d/%m/%Y"),
         year = lubridate::year(date))
ggplot(df) +
  geom_bar(aes(x = year, fill = Country_of_Residence)) +
  theme_bw() +
  scale_x_continuous(breaks = 2013:2024)

df2 <- geosentinel %>%
  select(source, Birth_Country, Country_of_Residence, Recent_Travel_Countries, Clinic_Visit_Date, Recent_Travel_Countries) %>%
  mutate(date = as.Date(Clinic_Visit_Date, format = "%d/%m/%Y"),
         year = lubridate::year(date)) %>%
  mutate(Country_of_Residence = ifelse(Country_of_Residence == "", "Unknown", Country_of_Residence)) %>%
  mutate(Country_of_Residence = ifelse(Country_of_Residence == "Portugal_Angola", "Portugal", Country_of_Residence))

df3 <- df2 %>%
  mutate(country_detection = case_when(source == "Libia" ~ "Portugal",
                                       Country_of_Residence != "Angola" ~ Country_of_Residence,
                                       Birth_Country != "Angola"~ Birth_Country,
                                       TRUE ~ "Germany"))

travel_time_plot <- ggplot(df3) +
  geom_bar(aes(x = year, fill = country_detection)) +
  theme_bw() +
  scale_x_continuous(breaks = 2013:2024)+
  labs(x = "Year", y = "DENV Positive Traveller Cases",
       fill = "Country of Residence") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

x <- df3 %>%
  filter(year(date) == "2013")
ggplot(x) +
  geom_bar(aes(x = date, fill = country_detection)) +
  theme_bw() +
  scale_x_date(limits = c(as.Date("2013-02-01", format = "%Y-%m-%d"),
                          as.Date("2013-07-31", format = "%Y-%m-%d")),
               breaks = "1 month", date_labels = "%b") +
  labs(x = "", y = "DENV Positive Traveller Cases",
       fill = "Country of Residence")

y <- df3 %>%
  filter(year(date) %in% c("2017", "2018", "2019"))
ggplot(y) +
  geom_bar(aes(x = date, fill = country_detection)) +
  theme_bw() +
  scale_x_date(limits = c(as.Date("2017-01-01", format = "%Y-%m-%d"),
                          as.Date("2019-12-31", format = "%Y-%m-%d"))) +
  labs(x = "Month of 2018", y = "DENV Positive Traveller Cases",
       fill = "Country of Residence") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


## Angola in and out by year flight volumes
total_out_angola_year <- angola %>%
  filter(orig_country == "Angola") %>%
  group_by(dest_country) %>%
  summarise(out_passengers = sum(passengers)) %>%
  dplyr::arrange(desc(out_passengers))

df4 <- df3 %>%
  group_by(country_detection) %>%
  summarise(sum = n()) %>%
  left_join(total_out_angola_year, by = c("country_detection" = "dest_country")) %>%
  filter(country_detection != "Unknown")

format_label <- function(x) {
  sapply(x, function(x) as.expression(bquote(10^.(x))))
}
average_travel_volume_vs_cases <- ggplot(df4, aes(x = log10(out_passengers / 10), y = sum, fill = country_detection)) +
  geom_point(size = 3, pch = 21) +
  theme_bw() +
  scale_x_continuous(
    labels = format_label,  # Use the custom function for the y-axis labels
    breaks = 2:5           # Set breaks at the actual y values
  ) +
  labs(x = "Average Yearly Air-Traffic Volume\nAngola -> Country of Residence",
       y = "Number of Detected Dengue\nCases in Travellers",
       fill = "Country of Residence")

geosentinel_plot <- cowplot::plot_grid(travel_time_plot + theme(legend.position = "none"), 
                                       average_travel_volume_vs_cases,
                                       ncol= 2,
                                       rel_widths = c(1, 1.3),
                                       align = "h", axis = "b", labels = c("B", "C"))

x <- cowplot::plot_grid(geosentinel_plot, flight_plot,
                   nrow = 2, align = "v", axis = "tblr")
ggsave(filename = "C:/Users/cw1716/Documents/Q_Drive_Copy2/Active_Research_Projects/angola_dengue/Fig5_FlightData_Figure.pdf", plot = x,
       width = 12, height = 9)


plot_legend <- cowplot::get_legend(average_travel_volume_vs_cases)

geosentinel_plot <- cowplot::plot_grid(travel_time_plot + theme(legend.position = "none"), 
                                       average_travel_volume_vs_cases + theme(legend.position = "none"),
                                       ncol= 2,
                                       rel_widths = c(1, 1.3),
                                       align = "h", axis = "b", labels = c("A", "B"))


x <- cowplot::plot_grid(geosentinel_plot, flight_plot,
                   nrow = 2, align = "v", axis = "tblr")


#######

total_out_angola_year <- angola %>%
  filter(orig_country == "Angola") %>%
  group_by(dest_country, year) %>%
  summarise(out_passengers = sum(passengers)) %>%
  dplyr::arrange(desc(out_passengers))

df3 <- df2 %>%
  group_by(year, Country_of_Residence) %>%
  summarise(sum = n()) %>%
  left_join(total_out_angola_year, by = c("Country_of_Residence" = "dest_country",
                                          "year")) %>%
  filter(Country_of_Residence != "Unknown") %>%
  group_by(Country_of_Residence) %>%
  filter(!is.na(out_passengers)) %>%
  summarise(years = n(),
            total_cases = sum(sum),
            total_passengers = sum(out_passengers),
            yearly_passengers = total_passengers / years)

better_average_travel_volume_vs_cases <- ggplot(df3, aes(x = log10(total_passengers), y = total_cases, fill = Country_of_Residence)) +
  geom_point(size = 3, pch = 21) +
  theme_bw() +
  scale_fill_discrete(name = "Country of\nResidence") +
  scale_x_continuous(
    labels = format_label,  # Use the custom function for the y-axis labels
    breaks = 2:6           # Set breaks at the actual y values
  ) +
  labs(x = "Total Air-Traffic Volume Angola -> Country of Residence\nIn Years With Detected Cases",
       y = "Total Number of Detected Dengue\nCases in Travellers")

cowplot::plot_grid(
  travel_time_plot + theme(legend.position = "none"), 
  better_average_travel_volume_vs_cases,
  ncol= 2,
  rel_widths = c(1, 1.4),
  align = "h", axis = "b", labels = c("A", "B"))

cor(log(df3$total_passengers), df3$total_cases)


plot(df3$sum, log(df3$out_passengers))

foi <- readRDS("angola_catalyticmodels_cleaned/foi_no_cross.rds")

df4 <- df2 %>%
  group_by(year) %>%
  summarise(sum = n()) %>%
  complete(year = 2013:2024, fill = list(sum = 0)) %>%
  left_join(foi, by = c("year"))

plot(df4$sum, df4$point_estimate)

table(df3$Country_of_Residence)
