# Resilience-of-spatial-structure-in-intestinal-organoids-revealed-by-spatial-transcriptomics
This is the repo of the code containing the necessary functions to reproduce the analyses and figure panels presented in our paper.


## Installation and Setup

1. Clone the repository.
2. Download the data from https://zenodo.org/records/17956270?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjczOWNlZDI0LTFlM2EtNGJjNy1iOTQwLTVlNDRjYTMzMjRlZCIsImRhdGEiOnt9LCJyYW5kb20iOiJkMGY2YjE2MjA4MDk0OWE2NTFiZjYxNzM4ZTI1YTgzZSJ9.Krydne1ABLYQe4DYJi13XDA6LZZ5mSfEpAeeTy32doKHx45VSEQEj8LUXBhq4CtHB1-Yfpmuih_lPA_jW-iaxg
3. Download from 'Moor, A. E., Harnik, Y., Ben-Moshe, S., Massasa, E. E., Rozenberg, M., Eilam, R., ... & Itzkovitz, S.
 (2018). Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis. Cell,
  175(4), 1156-1167.' 
  (a) table_A_LCM_TPM_values.tsv
  (b) table_D_zonation_reconstruction.tsv 
  and save them in an additional director 'in_vivo_villus_data'
  b. at this point the directories should be in the following structure:
  in a directory path of your choosing which will be defined in the variable HOME_DIR:
  ---in_vivo_villus_data 
    --table_A_LCM_TPM_values.tsv
    --table_D_zonation_reconstruction.tsv
  ---- sprinkled 
    -- 12hr
        --roi1
        --roi2
        --roi3
    -- 72hr
       --roi1
       --roi2
       --roi3
       --roi4
 ---- unperturbed
        --- monolayer_erosion
4. change the directory path for HOME_DIR in utils/constants.py to where you have saved these directories in the format explained above 
5. download packages using the requirements.txt file found in utils subdirectory
6. under paper/plotScripts exists the python file plotALL.py includes a function plot_all_figures() , call and run it in main.py 
7. the plots shown in the figures will be saved to paper/graphs/
to four different directions
    a. autonomous_zonation_figure_plots
    b. zonation_plasticity_plots
    c. neighborhood_zone_adoption_plots
    d. continuous_regenerative_figure_plots

"""
