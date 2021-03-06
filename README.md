# sparse-growth-curve
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) for the code.

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/) for other contents.

Analysis of bacterial growth curves with relatively low temporal resolution ([introduction and frequently asked questions](https://github.com/thrash-lab/sparse-growth-curve/blob/main/FAQ.md)). The code is written in Python in [**Google colab notebooks**](https://colab.research.google.com/notebooks/intro.ipynb) (compatible with [Jupyter notebook](https://jupyter.org/)).

1. To learn the **detailed method** of analyzing an individual growth curve (cell density vs. time), see [1_one_growth_curve_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb)
2. You can put multiple and related growth curves in one individual excel (.xlsx) file. As examples, you can see the ".xlsx" files in https://github.com/thrash-lab/sparse-growth-curve/tree/main/Growth_curve_data_example. For how an individual file is parsed, see [2_one_file_multiple_growth_curves_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb)
3. See [3_multiple_files_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb), for running multiple growth curve files in batch.


## Test run (To get started)
#### 1. Click the "Open In Colab" badge. In the screenshot below, click on (the red arrow) the file icon [1_growth_curve_analysis.ipynb](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb).
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/1_click_the_colab_badge.png" width="600">

#### 2. In the interface of the Colab environment, the first two things you can do are:
- Click the file icon (red arrow)- the file system will be shown. Here, we can ignore the "sample_data" folder. It is a folder created by Google Colab for their testing purposes.
- Hit the "Copy to drive" button (blue arrow)- a copy will then be made in your own Drive, and you can then make changes to the notebook if you so desire.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/2_file_folder_copy_to_drive.png" width="600">

- You can then find your notebook in the "Colab Notebooks" folder in your Google Drive.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/2_colab_notebook_folder.png" width="600">

#### 3. Select "Run all"- this will run all the code at once. All three notebooks are compatible with the "Run all" action.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/3_run_all.png" width="800">

- After running all the cells and hitting the refresh button (red arrow), you will see the output files. Then you can download the output files by right clicking the target (blue arrow).  
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/3_download_output.png" width="800">

#### 4. Besides running all the cells at once, you can also hit the play button to run each individual cell (e.g. for debugging purposes).
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/4_run_cell.png" width="800">

#### 5. In [2_one_file_multiple_growth_curves_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) and [3_multiple_files_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb), you need to upload your own Excel file(s) if you want to analyze your own data. Hit the upload file button and choose your data files.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/5_upload_files.png" width="800">


## Analyzing your own data
#### 1. For the format of your data, see the [**sample data**](https://github.com/thrash-lab/sparse-growth-curve/tree/main/Growth_curve_data_example). You **NEED** to have the **EXACT SAME** sheet tab names ("Data" and "Units"), column names in both sheets, and rows names in the "Units" sheet (MAKE SURE there are no extra spaces before and after those names). You **DO NOT NEED** to sort your data in a particular order. You **COULD** also add extra columns in the table. Those extra columns will not be read by the code.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/5_file_format.png" width="1200">

#### 2. Open [3_multiple_files_analysis.ipynb](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb)
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/6_open_multiple_files_colab.png" width="800">.

#### 3. You can delete the first code cell of the notebook. The first cell is for uploading the example data files through the command "wget". You do not need those example data files to analyze your own data.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/7_delete_first_cell.png" width="800">.

#### 4. Upload your own data files.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/8_upload_your_own_data_files.png" width="800">.

#### 5. Run all the code.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/9_run_all_multiple_files.png" width="800">.

#### 6. After executing all the code, there will be a .zip file created as output. Download this file- all the output tables and figures are inside.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/10_download_output.png" width="800">.

#### 7. If you get any errors, you can scroll down the notebook to "Run the script in batch" section and see which file is processed most recently. This is usually where an error will occur. Double check if the format of the file is correct or not.
<img src="https://github.com/thrash-lab/sparse-growth-curve/blob/main/image/11_check_error.png" width="800">.
