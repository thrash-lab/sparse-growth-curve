# sparse-growth-curve
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) for the code.

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/) for other contents.


Analysis of bacterial growth curves with relatively low temporal resolution ([introduction in detail](https://github.com/thrash-lab/sparse-growth-curve/blob/main/introduction.md)). The code is written in Python in [**Google colab notebooks**](https://colab.research.google.com/notebooks/intro.ipynb) (compatible with Jupyter notebook).

1. To learn how I parse an individual growth curve (cell density vs. time), see [1_one_growth_curve_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb)
2. You can put multiple and related growth curves in one individual excel (.xlsx) file. As examples, you can see the ".xlsx" files in https://github.com/thrash-lab/sparse-growth-curve/tree/main/Growth_curve_data_example. For how an individual file is parsed, see [2_one_file_multiple_growth_curves_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb)
3. See [3_multiple_files_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb), for running multiple growth curves files in batch.

## To get started
#### 1. Click the "Open In Colab" badge. In the screenshot below, I am clicking (the red arrow) the badge for the file [1_growth_curve_analysis.ipynb](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb).
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/1_click_the_colab_badge.png" width="600">

#### 2. In the interface of the Colab environment, the first two things you can do are:
- Click the file icon (red arrow), the file system will be shown. Here, we can ignore the "sample_data" folder. It is a folder created by Google Colab for their testing purposes.
- Hit the "Copy to drive" button (blue arrow), a copy will then be made in your own Drive, and you can then make changes to the notebook.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/2_file_folder_copy_to_drive.png" width="600">

- You can then find your notebook in the "Colab Notebooks" folder in your Google drive.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/2_colab_notebook_folder.png" width="600">

#### 3. Select "Run all", you can have all the code run at a time. And all three notebooks are compatible for "Run all" action.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/3_run_all.png" width="600">

- After run all, 
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/3_download_output.png" width="600">
