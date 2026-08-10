import os
import gzip
import numpy as np
import pandas as pd
import cv2
import json
import scanpy as sc


from optparse import OptionParser

parser=OptionParser(description='Generates Spatial h5ad')
parser.add_option('-i', '--input', '--in', default='.', help="BSTViewer_project/subdata/L13_heAuto/ [default: %default]")
parser.add_option('--png', help="BSTViewer_project/he_roi_small.png")
parser.add_option('--library', help="library id")
parser.add_option('--type', help="The options are 'S1000', 'S2000A', 'S2000B', 'S3000', 'cell_split'")
parser.add_option('-o', '--output', '--out', default='.', help="where to output the gene count matrix [default: %default]")
(opts, args)=parser.parse_args()

print(opts)

if not os.path.exists(opts.output):
        os.mkdir(opts.output)

adata = sc.read_10x_mtx(opts.input, var_names='gene_symbols',cache=True)
barcode_pos = pd.read_csv(f'{opts.input}/barcodes_pos.tsv.gz', compression='gzip', sep='\t', names=["barcode","array_col","array_row"],header=None)
barcode_pos['in_tissue']=1
barcode_pos=barcode_pos[['barcode','in_tissue','array_row','array_col']]
adata.obs=barcode_pos
adata.obs.index = adata.obs['barcode'].to_list()
adata.obs=adata.obs[['in_tissue','array_row','array_col']]

obsm=barcode_pos[['array_col','array_row']]
obsm= obsm.to_numpy()
adata.obsm["spatial"] = obsm

he_img = cv2.imread(opts.png)
he_img=he_img/255
library_id = opts.library
adata.uns["spatial"] = {library_id: {}}
adata.uns["spatial"][library_id]["images"] = {}
adata.uns['spatial'][library_id]['images']['hires'] = he_img.astype(np.float32)
adata.uns['spatial'][library_id]['use_quality']='hires'



def cal_zoom_rate(width, height, type):
    # 根据设备类型设定标准宽度
    if type in ["S1000", "S2000A", "S2000B", "S3000"]:
        std_width = 1000
    elif type == "cell_split":
        std_width = 20000
    else:
        raise ValueError(f"Unsupported device type: {type}")

    # 计算标准高度
    if type == "S3000":
        std_height = std_width / (42 * 46) * (43 * 52 * np.sqrt(3) / 2.0)
    elif type == "S2000A":
        std_height = std_width / (76 * 31) * (75 * 36 * np.sqrt(3) / 2.0)
    elif type == "S2000B":
        std_height = std_width / (101 * 31) * (134 * 36 * np.sqrt(3) / 2.0)
    else:  # 处理S1000/cell_split及其他默认情况
        std_height = std_width / (46 * 31) * (46 * 36 * np.sqrt(3) / 2.0)

    # 计算缩放比例
    aspect_ratio = std_width / std_height
    target_ratio = width / height
    
    if aspect_ratio > target_ratio:
        scale = width / std_width
    else:
        scale = height / std_height
    
    return scale

zoom_scale = cal_zoom_rate(he_img.shape[1], he_img.shape[0], opts.type)
adata.uns['spatial'][library_id]['scalefactors']= {"spot_diameter_fullres": zoom_scale, "tissue_hires_scalef": zoom_scale,
                                                   "fiducial_diameter_fullres": zoom_scale, "tissue_lowres_scalef": zoom_scale}
sc.pl.spatial(adata, img_key="hires")

adata.write(f'{opts.output}/adata.h5ad')

