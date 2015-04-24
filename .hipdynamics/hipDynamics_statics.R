
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+   _______ __         ______                             __               +")
print("+  |   |   |__|.-----.|      |--.--.-----.---.-.--------.|__|.----.-----.  +")
print("+  |       |  ||  _  ||  --  |  |  |     |  _  |        ||  ||  __|__ --|  +")
print("+  |___|___|__||   __||_____/|___  |__|__|___._|__|__|__||__||____|_____|  +")
print("+              |__|          |_____|                                       +")
print("+                                                                          +")
print("+ TITLE:   HipDynamics - An analysis to deduce cell population dynamics    +")
print("+ VERSION: 1.2                                                             +")
print("+ AUTHOR:  Maximilian Kerz (kerz.maximilian@gmail.com)                     +")
print("+                                                                          +")
print("+ ACKNOWLEDGEMENTS: Amos Folarin (amosfolarin@gmail.com)                   +")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")                                  

# ------------------------------------ Statics -------------------------------------

# index
img_no_idx <- "ImageNumber"
file_idx <- "Image_FileName_Phase"
well_idx <- "Image_Metadata_Well"
order_idx <- "Order_ID"
fn_idx <- "FN_conc"
line_idx <- "Cell_Line"
obj_img_no_idx <- "ImageNumber"

img_count_cellsFinal_idx <- "Image_Count_CellsFinal"
img_count_nuc1_idx <- "Image_Count_Nuclei1"
img_count_nuc2_idx <- "Image_Count_Nuclei2"
img_count_nuc3_idx <- "Image_Count_Nuclei3"
img_count_nuc4_idx <- "Image_Count_Nuclei4"
img_cellsFinal_area_idx <- "Mean_CellsFinal_AreaShape_Area"
img_cellsFinal_compact_idx <- "Mean_CellsFinal_AreaShape_Compactness"
img_cellsFinal_accent_idx <- "Mean_CellsFinal_AreaShape_Eccentricity"
img_cellsFinal_maxRad_idx <- "Mean_CellsFinal_AreaShape_MaximumRadius "
img_cellsFinal_meanRad_idx <- "Mean_CellsFinal_AreaShape_MeanRadius"
img_cellsFinal_medRad_idx <- "Mean_CellsFinal_AreaShape_MedianRadius"

obj_area_idx <- "CellsFinal_AreaShape_Area"
obj_compact_idx <- "CellsFinal_AreaShape_Compactness"
obj_eccent_idx <- "CellsFinal_AreaShape_Eccentricity"
obj_euler_idx <- "CellsFinal_AreaShape_EulerNumber"
obj_formFactor_idx <- "CellsFinal_AreaShape_FormFactor"
obj_majorAxisLength_idx <- "CellsFinal_AreaShape_MajorAxisLength"
obj_maxRad_idx <- "CellsFinal_AreaShape_MaximumRadius"
obj_meanRad_idx <- "CellsFinal_AreaShape_MeanRadius"
obj_medRad_idx <- "CellsFinal_AreaShape_MedianRadius"
obj_minAxisLength_idx <- "CellsFinal_AreaShape_MinorAxisLength"
obj_perimeter_idx <- "CellsFinal_AreaShape_Perimeter"


obj_index_list <- c(obj_area_idx
                    #obj_compact_idx,
                    #obj_eccent_idx,
                    #obj_euler_idx,
                    #obj_formFactor_idx,
                    #obj_majorAxisLength_idx,
                    #obj_maxRad_idx,
                    #obj_meanRad_idx,
                    #obj_medRad_idx,
                    #obj_minAxisLength_idx,
                    #obj_perimeter_idx
)

thresh_last <- threshold[length(threshold)]
fields_list <- length(obj_index_list)

# Experiment's-id-storage-variable (from plate results)
experiments <- c()
cell_lines <- c()
# number of hours of live-imaging
img_hours <- 24
img_hours_suggest=c(24,23,22,21,20)
img_hours_lowest <- 20

# Prompt
aPrompt <- "[HipDynamics]"

