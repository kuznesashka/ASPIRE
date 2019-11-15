# ASPIRE

Automated interictal spike detection and source localization in magnetoencephalography using independent components analysis and spatio-temporal clustering

## Requirements 
* Brainstorm 3
* FieldTrip 


## Folders structue

```bash
ROOT
├─── ASPIRE
│    ├── detections
│    │   ├── ICA_detections_B1C2_ii_run1_raw_tsss_mc_grad.mat
│    │   ├── ICA_detections_B1C2_ii_run1_raw_tsss_mc_mag.mat
│    │   ├── Manual_spikes_B1C2_ii_run1_raw_tsss_mc.csv
│    │   ├── Templates_B1C2_ii_run1_raw_tsss_mc_grad.csv
│    │   └── Templates_B1C2_ii_run1_raw_tsss_mc_mag.csv
│    ├── plots
│    │   ├── B1C2.bmp
│    │   ├── clusters
│    │   │   ├── Cluster_ICA_based_grad_1.bmp
│    │   │   ├── Cluster_ICA_based_mag_23.bmp
│    │   │   ├── Cluster_SpyCir_based_mag_1.bmp
│    │   │   └── Cluster_SpyCir_based_mag_30.bmp
│    │   ├── overlap_grad.png
│    │   ├── overlap_grad_ASPIRE_1_SPC_18.png
│    │   ├── overlap_grad_ASPIRE_1_SPC_60.png
│    │   ├── overlap_grad_ASPIRE_2_SPC_18.png
│    │   ├── overlap_grad_ASPIRE_7_SPC_29.png
│    │   ├── overlap_mag.png
│    │   ├── overlap_mag_ASPIRE_1_SPC_40.png
│    │   ├── overlap_mag_ASPIRE_2_SPC_37.png
│    │   └── overlap_mag_ASPIRE_6_SPC_16.png
│    └── results
│        ├── cluster_out_ICA_based_grad.csv
│        ├── cluster_out_ICA_based_mag.csv
│        ├── cluster_out_SpyCir_based_grad.csv
│        ├── cluster_out_SpyCir_based_mag.csv
│        ├── cluster_out_overlap_grad.csv
│        ├── cluster_out_overlap_grad.mat
│        ├── cluster_out_overlap_mag.csv
│        ├── cluster_out_overlap_mag.mat
│        ├── cluster_out_time_only_ICA_based_grad.mat
│        ├── cluster_out_time_only_ICA_based_mag.mat
│        ├── cluster_out_time_only_SpyCir_based_grad.mat
│        ├── cluster_out_time_only_SpyCir_based_mag.mat
│        ├── cluster_out_time_only_visual_grad.mat
│        ├── cluster_out_time_only_visual_mag.mat
│        ├── cluster_out_visual_grad.csv
│        ├── cluster_out_visual_mag.csv
│        ├── results_ICA_based_grad.mat
│        ├── results_ICA_based_mag.mat
│        ├── results_SpyCir_based_grad.mat
│        ├── results_SpyCir_based_mag.mat
│        ├── results_visual_grad.mat
│        ├── results_visual_mag.mat
│        ├── sources_ICA_based_grad.mat
│        ├── sources_ICA_based_mag.mat
│        ├── sources_SpyCir_based_grad.mat
│        ├── sources_SpyCir_based_mag.mat
│        ├── sources_visual_grad.mat
│        └── sources_visual_mag.mat
│            
└── anat_and_data_folder    
	├── anat
	│   └── B1C2
	│       ├── brainstormsubject.mat
	│       ├── subjectimage_T1.mat
	│       ├── tess_aseg.mat
	│       ├── tess_cortex_mid_high.mat
	│       ├── tess_cortex_mid_low.mat
	│       ├── tess_cortex_pial_high.mat
	│       ├── tess_cortex_pial_low.mat
	│       ├── tess_cortex_pialcereb_low.mat
	│       ├── tess_cortex_white_high.mat
	│       ├── tess_cortex_white_low.mat
	│       └── tess_head_mask.mat
	└── data
	    ├── @default_study
	    │   └── brainstormstudy.mat
	    ├── @inter
	    │   └── brainstormstudy.mat
	    └── B1C2
	        ├── @default_study
	        │   ├── brainstormstudy.mat
	        │   ├── channel_vectorview306_acc1.mat
	        │   └── headmodel_surf_os_meg.mat
	        └── B1C2_ii_run1_raw_tsss_mc_art_corr
	            ├── brainstormstudy.mat
	            ├── channel_vectorview306_acc1.mat
	            ├── data_block001.mat
	            ├── data_block002.mat
	            ├── data_block003.mat
	            └── headmodel_surf_os_meg.mat
	
```


## Literature
* Ossadtchi et al., 2004 (https://www.sciencedirect.com/science/article/abs/pii/S138824570300405X)
