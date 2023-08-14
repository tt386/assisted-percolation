#!/bin/bash

path_dst=$1


Sources=(
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_1B_i_Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.300/Background.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_1B_ii_Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.600/Background.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_1C_Height_1500_Width_1000_Fitness_0.950_MutProb_0.00100_Radius_50_Random_AreaFraction_0.500/AField_00075.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_1C_Height_1500_Width_1000_Fitness_0.950_MutProb_0.00100_Radius_50_Random_AreaFraction_0.500/AField_00150.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_1C_Height_1500_Width_1000_Fitness_0.950_MutProb_0.00100_Radius_50_Random_AreaFraction_0.500/SystemEvolution.mp4

    Scripts_OneNode/ClusterHeightStatistics/SaveFiles/SampleOutput/Rough_Vs_Flat.pdf
    Scripts_HPC/InfectionProbability/SaveFiles/Fig_Height_280_Width_180_Repeats_10000_Radius_60_BottomBuffer_100_FlatInfection_RoughFront/InfectionProb.pdf
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00011.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00019.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00025.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/SystemEvolution.mp4
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00011.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00019.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00025.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/SystemEvolution.mp4

    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00025.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00042.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00063.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00025.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00042.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/AField_00063.png
    Sketch/SymmetricEscapeRegion/SymmetricPatch.pdf
    Scripts_HPC/EscapeRegion_HeightVsFitness/SaveFiles/Fig_Height_1000_Width_150_Repeats_1000_Radius_50_Min_Fitness_0.700_Max_Fitness_0.990_num_Fitness_100_FlatInfection_RoughFront_Random_Radius_50/Median_CapHeight.pdf
    Scripts_HPC/EscapeRegion_HeightVsAngle/SaveFiles/Fig_Height_1000_Width_150_Repeats_1000_Radius_50_Fitness_0.950_FlatInfection_RoughFront_Random_Radius_50/Median_CapHeight.pdf

    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/SystemEvolution.mp4
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/SystemEvolution.mp4
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.700/AField_00210.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.700/SystemEvolution.mp4
    Scripts_HPC/PhaseDiagram_Random/SaveFiles/SampleOutput/pcolor_Domination.eps

    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/SystemEvolution.mp4
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/SystemEvolution.mp4
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_6.000_AreaFraction_0.800/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_6.000_AreaFraction_0.800/SystemEvolution.mp4
    Scripts_HPC/PhaseDiagram_Lattice/SaveFiles/SampleOutput/pcolor_Domination.eps

    Scripts_OneNode/SingleSystemEvolution/SaveFiles/ZoomedRoughFront_Height_30_Width_30_Fitness_0.800_MutProb_0.01000_Radius_5_SingleObstacle/AField_00003.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.02000_Radius_40_Random_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/AField_00200.png
    Scripts_OneNode/SingleSystemEvolution/SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.700/AField_00210.png

    Scripts_OneNode/FlatFrontEvolution/SaveFiles/SampleOutput/AField_00003.png
    Scripts_OneNode/ClusterHeightStatistics/SaveFiles/SampleOutput/WidthGivenHeightPlots/Height_0010.eps
    Scripts_OneNode/ClusterHeightStatistics/SaveFiles/SampleOutput/WidthGivenHeightPlots/Height_0020.eps
    Scripts_OneNode/ClusterHeightStatistics/SaveFiles/SampleOutput/WidthGivenHeightPlots/Height_0030.eps

    Sketch/EscapeRegionCalculation/EscapeCalculation.pdf
    Sketch/EscapeRegionCalculation_Over/EscapeCalculation.pdf

    Sketch/Lattice_Orientation/Horizontal.eps
    Sketch/Lattice_Orientation/Vertical.eps
    Sketch/Lattice_ShortTerm/ShortTerm.pdf
    Sketch/Lattice_LongTerm/LongTerm.pdf
)


Names=(
    Fig_1_B_i.png
    Fig_1_B_ii.png
    Fig_1_C_i.png
    Fig_1_C_ii.png
    Vid_1_C.mp4

    Fig_2_A.pdf
    Fig_2_B.pdf
    Fig_2_C_i_a.png
    Fig_2_C_i_b.png
    Fig_2_C_i_c.png
    Vid_2_C_i.mp4
    Fig_2_C_ii_a.png
    Fig_2_C_ii_b.png
    Fig_2_C_ii_c.png
    Vid_2_C_ii.mp4   

    Fig_3_A_i_a.png
    Fig_3_A_i_b.png
    Fig_3_A_i_c.png
    Fig_3_A_ii_a.png
    Fig_3_A_ii_b.png
    Fig_3_A_ii_c.png
    Fig_3_B.pdf
    Fig_3_C.pdf
    Fig_3_D.pdf

    Fig_4_A_i.png
    Vid_4_A_i.mp4
    Fig_4_A_ii.png
    Vid_4_A_ii.mp4
    Fig_4_A_iii.png   
    Vid_4_A_iii.mp4
    Fig_4_B.eps
    
    Fig_5_A_i.png
    Vid_5_A_i.mp4
    Fig_5_A_ii.png
    Vid_5_A_ii.mp4
    Fig_5_A_iii.png
    Vid_5_A_iii.mp4
    Fig_5_B.eps

    SuppFig_1_B.png
    SuppFig_1_C_i.png
    SuppFig_1_C_ii.png
    SuppFig_1_C_iii.png
    SuppFig_1_C_iiii.png

    SuppFig_2_A.png
    SuppFig_2_B_i.eps
    SuppFig_2_B_ii.eps
    SuppFig_2_B_iii.eps

    SuppFig_3_A.pdf
    SuppFig_3_B.pdf

    SuppFig_5_Ai.eps
    SuppFig_5_Aii.eps
    SuppFig_5_B.pdf
    SuppFig_5_C.pdf
)





for i in "${!Sources[@]}"; do
    basename "${Sources[$i]}"
    f="${Names[$i]}"
    echo $filename
    file_dst="${path_dst}/${f}"
    
    echo $file_dst

    cp "${Sources[$i]}" "$file_dst"
    echo cp "${Sources[$i]}" "$file_dst"
done
