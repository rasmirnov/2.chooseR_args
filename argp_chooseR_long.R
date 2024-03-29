# PATH to INPUT: "/mnt/tank/scratch/rasmirnov/code/Time_benchmarking/Datasets/processed/obj/obj_tutorial.RData"
# PATH to OUTPUT: "/mnt/tank/scratch/rasmirnov/code/chooseR/results/pbmc3k/argparse_res3"

source("/mnt/tank/scratch/rasmirnov/code/chooseR/scripts/pipeline/pipeline.R") 
library(Seurat)
library(microbenchmark)
library(ggplot2)
library(argparse)
library(R6)
library(findpython)
library(jsonlite)
`%>%` <- magrittr::`%>%`

parser <- ArgumentParser(description='Time benchmarking of the chooseR tool')
parser$add_argument('-d', 
                    '--data',
                    type = 'character',
                    help = 'path to rdata object')
parser$add_argument('-r', 
                    '--res',
                    type = 'character',
                    help = 'path to output directory')
parser$add_argument('-i',                                   # immune.combined.sct
                    '--obj_integ',
                    type = 'character',
                    help = 'name of a seurat_object')
parser$add_argument('-o',                                   # ifnb, stestis, etc
                    '--obj_name',
                    type = 'character',
                    help = 'name of a seurat_object')


argss <- parser$parse_args()                          # ?separation of args inside the container
argss
load(argss$data)
object <- get(argss$obj_integ)
obj_name <- argss$obj_name
results_path <- argss$res
# results_path <- '/mnt/tank/scratch/rasmirnov/code/chooseR/results/stestis'
# resolutions <- 0.05
resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 5)                     # 0.05 was removed

setwd(results_path)

# Run pipeline
start.time <- Sys.time()

chooseR <- function(seurat_obj, 
                    resolutions,
                    n, 
                    results_path, 
                    npcs = 50,
                    assay = 'SCT', 
                    reduction = 'pca') {
  results_path <- sprintf("%s/n_%s", results_path, n)
  message(sprintf("n = %s, results directory: %s", n, results_path))
  for (res in resolutions) {
    message(paste0("Clustering ", res, "..."))
    message("\tFinding ground truth...")
    
    clusters <- seurat_obj[[glue::glue("integrated_snn_res.{res}")]]
    
    # Now perform iterative, sub-sampled clusters
    results <- multiple_cluster(
      seurat_obj,
      n = n,
      size = 0.8,
      npcs = npcs,
      res = res,
      reduction = reduction,
      assay = assay
    )
    
    # Now calculate the co-clustering frequencies
    message(paste0("Tallying ", res, "..."))
    columns <- colnames(dplyr::select(results, -cell))
    mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
    i <- 1 # Counter
    for (col in columns) {
      message(paste0("\tRound ", i, "..."))
      mtchs <- Reduce("+", list(
        mtchs,
        find_matches(col, df = results)
      ))
      i <- i + 1
    }
    
    message(paste0("Scoring ", res, "..."))
    mtchs <- dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
    )
    
    # Now calculate silhouette scores
    message(paste0("Silhouette ", res, "..."))
    sil <- cluster::silhouette(
      x = as.numeric(as.character(unlist(clusters))),
      dmatrix = (1 - as.matrix(mtchs))
    )
    if (!file.exists(paste0(results_path, "silhouette_", res, ".rds"))) {
      saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
    }
    # Finally, calculate grouped metrics
    message(paste0("Grouping ", res, "..."))
    grp <- group_scores(mtchs, unlist(clusters))
    if (!file.exists(paste0(results_path, "frequency_grouped_", res, ".rds"))) {
      saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
    }
    
    sil <- group_sil(sil, res)
    if (!file.exists(paste0(results_path, "silhouette_grouped_", res, ".rds"))) {
      saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
    }
  }
  if (!file.exists(paste0(results_path, "clustered_data.rds"))) {
    saveRDS({{seurat_obj}}, paste0(results_path, "clustered_data.rds"))
  }
}


# Microbenchmarking + Function call
mbm <- microbenchmark(
  pbmc3k_n100 = chooseR(object, resolutions, n=100, results_path = results_path, npcs = 50),        # n=100
  times = 1
)

#### Extra timing
end.time <- Sys.time() - start.time

# Create silhouette plot
# Read in scores and calculate CIs
# set preferable 'n' which you want research:  
n=100                            #n=100
results_path <- sprintf("%s/n_%s", results_path, n)                #! очень важный элемент
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()
meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

##! all data except 0.05 resolution:
# meds <- filter(slice(meds, 2:10))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
    dplyr::filter(med >= threshold) %>%
    dplyr::arrange(n_clusters) %>%
    tail(n = 1) %>%
    dplyr::pull(res)
)

# And plot!
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  geom_hline(aes(yintercept = threshold), colour = "blue") +
  geom_vline(aes(xintercept = choice), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = "silhouette_distribution.png",
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Finally, a dot plot of silhouette scores to help identify less robust clusters
# The initial pipe is to order the clusters by silhouette score
scores %>%
  dplyr::filter(res == choice) %>%
  dplyr::arrange(dplyr::desc(avg_sil)) %>%
  dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
  ggplot(aes(factor(cluster), avg_sil)) +
  geom_point() +
  scale_x_discrete("Cluster") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_grid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0("silhouette_point_plot_", choice, ".png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# # save.image('pbmc3k_chooseR.RData')
save.image(file = paste0(obj_name, '_chooseR.RData'))
