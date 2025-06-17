# overlap_disgenet_genecards_opentargets

# This project is to create an overlap between three databases which are DisGeNET, Gene Cards and OpenTargets to find shared genes between Ankylosing Spondylitis (AS) and Psoriasis (PsO)

# Step 1: Data downloaded in form of TSV for AS and PsO separately in DisGeNET, Gene Cards and OpenTargets

DisGeNET --> https://disgenet.com/ (Limited data for unsubscriber, may need to search manually WITHOUT login until certain threshold you want)
Gene Cards --> https://www.genecards.org/ (Free)
OpenTargets --> https://platform.opentargets.org/ (Free)

Ensure not duplicates in every datasets

# Step 2: Overlap datasets in R

Change the paths of the datasets to your own computer's paths that contain downloaded datasets

If haven't download R, can download R and Rstudio at https://posit.co/download/rstudio-desktop/

## Load the code and generate Venn Diagram for the overlap between DisGeNET, Gene Cards and OpenTargets
