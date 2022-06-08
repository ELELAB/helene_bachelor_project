# Improving the workflows for analyses of gene expression data from The Cancer Genome Atlas

The program files (R scripts - files with extension .R) in this repository
are licensed under the terms of the GNU General Public License (see LICENSE
file). You can redistribute them and/or modify them under the terms of the
GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
These programs are distributed in the hope that they will be useful,but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License
along with this program. If not, see https://www.gnu.org/licenses/.

Each folder has a README file to explain the purpose of each program and
output files. The directory organization are inspired by the Cancer Systems
Biology groups server organization.

This thesis uses data from 31 cancer types originating from the cancer genome
atlas (TCGA), located in the /data. There are three main parts of the analysis:

All the cancer data are used to investigate an appropriate method for
handling sample replicates located in /Expression/cancer_type/Technical_replicates.

The breast cancer (BRCA) data set are used to verify a new generated
normalization table located. The new table are located in
/Expression/get_GLGC_content while the BRCA analysis are performed in
/Expression/cancer_type/BRCA_getGLGC and /Expression/cancer_type/DEA_getGLGC

A case study on Uterine Corpus Endometrial Carcinoma (UCEC) are performed
with the new pipeline including handling replicates and normalization by the
new table in the folder /Expression/cancer_type/Casestudy. Followed up by
differential expression, enrichment analysis and visualizations for both UCEC
and subtype specific in the subfolders to /Casestudy.

