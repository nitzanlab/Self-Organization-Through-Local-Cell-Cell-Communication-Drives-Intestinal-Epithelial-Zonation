"""
This file includes all the imports needed to run the code.
In order to install them on your machine, or install to your conda environment
run in the terminal : pip install -r requirements.txt
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
from scipy import stats
import numpy as np
from scipy.stats import entropy
import novosparc
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr

import glob
from scipy.signal import convolve2d
from scipy.spatial import cKDTree
from sklearn.neighbors import KDTree
from scipy.stats import ks_2samp
from sklearn.metrics.pairwise import cosine_similarity
from matplotlib.lines import Line2D
import pickle
import matplotlib.cm as cm

from skimage import morphology
from skimage.morphology import binary_closing, disk, binary_opening, erosion


from scipy.ndimage import gaussian_filter
from skimage import morphology
from skimage.morphology import binary_opening, disk
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from sklearn.neighbors import BallTree

from skimage.measure import regionprops, label