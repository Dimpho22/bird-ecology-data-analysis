# Bird Niche Breadth and Niche Overlap Analysis
# Data: bird counts recorded at 28 Wits University habitat patches.
# Important: the habitat patch names and counts below come from the supplied dataset.
# The script treats each sampled patch as equally available when calculating
# Smith's niche-breadth index because patch area/availability data were not supplied.
# Interpret niche breadth cautiously: equal availability is an assumption, not measured data.

library(tidyverse)
library(janitor)

bird_data <- read_csv(file.choose(), show_col_types = FALSE) |>
  clean_names()

species_cols <- setdiff(names(bird_data), "habitat_patch")

stopifnot(nrow(bird_data) == 28)
stopifnot(length(species_cols) == 6)
bird_data[is.na(bird_data)] <- 0

bird_long <- bird_data |>
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "species",
    values_to = "count"
  ) |>
  mutate(
    species = str_replace_all(species, "_", " "),
    species = str_to_title(species)
  )

get_mode <- function(x) {
  counts <- table(x)
  as.numeric(names(counts)[counts == max(counts)])
}

species_summary <- bird_long |>
  group_by(species) |>
  summarise(
    total_count = sum(count),
    mean_per_patch = mean(count),
    median_per_patch = median(count),
    mode_per_patch = paste(get_mode(count), collapse = ", "),
    standard_deviation = sd(count),
    standard_error = sd(count) / sqrt(n()),
    patches_present = sum(count > 0),
    patches_absent = sum(count == 0),
    .groups = "drop"
  ) |>
  arrange(desc(total_count))

print(species_summary)
write_csv(species_summary, "species_summary.csv")

patch_summary <- bird_long |>
  group_by(habitat_patch) |>
  summarise(total_birds = sum(count), .groups = "drop") |>
  arrange(desc(total_birds))

print(patch_summary)
write_csv(patch_summary, "habitat_patch_summary.csv")

histogram_plot <- ggplot(bird_long, aes(x = count)) +
  geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = "Bird counts across habitat patches",
    x = "Count recorded in a habitat patch",
    y = "Number of patches"
  ) +
  theme_minimal()

print(histogram_plot)
ggsave("bird_count_histograms.png", histogram_plot, width = 10, height = 7, dpi = 300)

boxplot_plot <- ggplot(bird_long, aes(x = species, y = count)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    title = "Distribution of bird counts by species",
    x = "Species",
    y = "Count per habitat patch"
  ) +
  theme_minimal()

print(boxplot_plot)
ggsave("bird_count_boxplots.png", boxplot_plot, width = 9, height = 6, dpi = 300)

species_totals <- colSums(bird_data[species_cols], na.rm = TRUE)
print(species_totals)

chi_square_result <- chisq.test(species_totals, p = rep(1 / length(species_totals),
                                                        length(species_totals)))
print(chi_square_result)

chi_square_output <- tibble(
  test = "Chi-square goodness-of-fit: equal observed counts across species",
  chi_square = unname(chi_square_result$statistic),
  degrees_of_freedom = unname(chi_square_result$parameter),
  p_value = chi_square_result$p.value
)
write_csv(chi_square_output, "chi_square_result.csv")

smith_breadth <- function(counts, availability = rep(1 / length(counts), length(counts))) {
  if (sum(counts) == 0) return(NA_real_)
  p <- counts / sum(counts)
  sum(sqrt(p * availability))
}

set.seed(2442463)
n_boot <- 2000

breadth_results <- map_dfr(species_cols, function(sp) {
  counts <- bird_data[[sp]]
  estimate <- smith_breadth(counts)
  if (sum(counts) == 0) {
    return(tibble(species = sp, smith_b = NA_real_,
                  ci_2_5 = NA_real_, ci_97_5 = NA_real_))
  }

  boot_values <- replicate(n_boot, {
    idx = sample(seq_along(counts), replace = TRUE)
    sampled_counts = counts[idx]
    smith_breadth(sampled_counts, rep(1 / length(sampled_counts), length(sampled_counts)))
  })

  tibble(
    species = sp,
    smith_b = estimate,
    ci_2_5 = unname(quantile(boot_values, 0.025, na.rm = TRUE)),
    ci_97_5 = unname(quantile(boot_values, 0.975, na.rm = TRUE))
    )
}) |>
  mutate(species = str_to_title(str_replace_all(species, "_", " ")))

print(breadth_results)
write_csv(breadth_results, "niche_breadth_smith_index.csv")

pianka_overlap <- function(x, y) {
  if (sum(x) == 0 || sum(y) == 0) return(NA_real_)
  p <- x / sum(x)
  q <- y / sum(y)
  sum(p * q) / sqrt(sum(p^2) * sum(q^2))
}

species_pairs <- combn(species_cols, 2, simplify = FALSE)

overlap_results <- map_dfr(species_pairs, function(pair) {
  x <- bird_data[[pair[1]]]
  y <- bird_data[[pair[2]]]
  tibble(
    species_1 = str_to_title(str_replace_all(pair[1], "_", " ")),
    species_2 = str_to_title(str_replace_all(pair[2], "_", " ")),
    pianka_overlap = pianka_overlap(x, y)
  )
}) |>
  arrange(desc(pianka_overlap))

print(overlap_results)
write_csv(overlap_results, "niche_overlap_pianka.csv")

overlap_plot <- ggplot(overlap_results,
                       aes(x = reorder(paste(species_1, "vs", species_2),
                                       pianka_overlap),
                           y = pianka_overlap)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Pairwise niche overlap (Pianka's index)",
    x = "Species pair",
    y = "Pianka overlap (0–1)"
  ) +
  theme_minimal()

print(overlap_plot)
ggsave("pianka_niche_overlap.png", overlap_plot, width = 9, height = 7, dpi = 300)

cat("\n--- QUICK INTERPRETATION ---\n")
cat("Most frequently recorded species by total count:\n")
print(species_summary |> slice_max(total_count, n = 1, with_ties = TRUE))
cat("\nLeast frequently recorded species by total count:\n")
print(species_summary |> slice_min(total_count, n = 1, with_ties = TRUE))
cat("\nHabitat patch(es) with the highest total count:\n")
print(patch_summary |> slice_max(total_birds, n = 1, with_ties = TRUE))
cat("\nHabitat patch(es) with the lowest total count:\n")
print(patch_summary |> slice_min(total_birds, n = 1, with_ties = TRUE))
cat("\nHighest Pianka overlap:\n")
print(overlap_results |> slice_max(pianka_overlap, n = 1, with_ties = TRUE))
cat("\nNote: these results describe the supplied observations only.\n")
cat("The dataset contains many zero counts; sampling effort and detection may vary.\n")
