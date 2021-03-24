# 1. Median Country Ages -------------------------------------------------------


# url <- "https://en.wikipedia.org/wiki/List_of_countries_by_median_age"
# ages <- url %>% read_html() %>%
#   html_node(css = "body #content #bodyContent #mw-content-text .mw-parser-output table") %>%
#   html_table(fill = TRUE) %>%
#   tail(-1) %>%
#   select(c(1,3)) %>%
#   setNames(c("country","age")) %>%
#   mutate(iso3c = countrycode::countrycode(country, "country.name.en", "iso3c",
#                                           custom_match = c("Kosovo" = "KSV",
#                                                            "Virgin Islands" = "VIR",
#                                                            "Saint Martin" = "MAF")))

url <- "https://www.worldometers.info/world-population/population-by-country/"
iso3c_ages <- read_html(url) %>%
  html_table() %>%
  .[[1]] %>%
  select(c("Country (or dependency)", "Med. Age")) %>%
  setNames(c("country", "age")) %>%
  mutate(iso3c = countrycode::countrycode(
    country, "country.name.en", "iso3c",
    custom_match = c("Kosovo" = "KSV",
                     "Virgin Islands" = "VIR",
                     "Saint Martin" = "MAF",
                     "Channel Islands" = "CHI",
                     "Micronesia" = "FSM"))) %>%
  filter(iso3c %in% squire::population$iso3c) %>%
  mutate(age = as.numeric(age))
usethis::use_data(iso3c_ages, overwrite = TRUE)

# End --------------------------------------------------------------------------
