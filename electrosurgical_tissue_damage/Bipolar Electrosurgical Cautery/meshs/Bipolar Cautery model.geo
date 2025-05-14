SetFactory("OpenCASCADE");
Merge "Bipolar Cautery model.step";
//+
Physical Surface("electrode_1_border", 85) = {8, 1, 6, 3, 5, 4, 7};
//+
Physical Surface("electrode_2e_border", 86) = {11, 16, 10, 9, 15, 14, 12, 13};
//+
Physical Surface("electrode_1_border", 85) += {2};
//+
Physical Surface("Free_convection", 87) = {19, 22, 18, 20, 17};
//+
Physical Surface("T_0", 88) = {21, 24};
//+
Physical Surface("Normal", 89) = {23};
//+
Physical Volume("Electrode_1", 90) = {1};
//+
Physical Volume("Electrode_2", 91) = {2};
//+
Physical Volume("Tissue", 92) = {3};
//+
Physical Surface("electrode_1_tissue", 93) = {28, 27, 26, 25};
//+
Physical Surface("electrode_2_tissue", 94) = {32, 31, 30, 33};
//+
Physical Surface("electrode_1_tissue", 93) += {29};
//+
Physical Surface("electrode_2_tissue", 94) += {34};
