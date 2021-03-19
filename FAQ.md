#### Q1: What is sparse-growth-curve?
- A1: The sparse-growth-curve package is for analyzing multiple growth curves (can be under low temporal resolution) under multiple conditions and multiple strains with multiple replicates. You can refer to [sparse-growth-curve: A Computation Pipeline for Parsing Cellular Growth Curves under Low Temporal Resolution](https://docs.google.com/document/d/11uvXfO1qYBQnJ3qEPmVq4zH_UmXkScBMF8pW1s1PJgM/edit?usp=sharing) to get a brief idea of why and how wer disigned sparse-growth-curve package.

#### Q2: What is your method in detail to parse and characterize a single growth curve?
- A2: When you have a growth curve, namely an array of cell densities $$X_i$$ and the corresponding aquisition time $$t_i$$


#### Q3: Why are you using both [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) and [![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)?
- A3: We fully embrace the open-source software initiative, and therefore, **the code** are all under the **MIT license**. On the other hand, we want to have our code being versatile and can be run and understood naively. Therefore, we have written a lot of instructions of both our methods and the code. Those content are under the CC-By-NC license.

#### Q4: Why are you using [Google Colab notebook](https://colab.research.google.com/notebooks/intro.ipynb)?
- A4: On one hand, we want our methods to be run smoothly by users who do not have any programming background. Google Colab is a programming interface that can be processed completely online without any extra configurations. On the other hand, we also want our code to be versatile. Therefore instead of having our code encapsulated as black boxes of functions, we would rather have it presented step by step and layer by layer. In this case, we believe Google Colab is a good format for communication.

#### Q5: Why do you have 3 different notebooks? What are each of them for?
- A5: The sparse-growth-curve is designed for reading experimental data files, and parse the growth curves within them. Therefore, there are multiple layers of iterations involved:
    - (1) First, we need to iterate through every file in the folder.
    - (2) For each file in the folder, we need to parse each individual growth curve out.
    - (3) For each individual growth curve, we need to apply the sparse-growth-curve method and characterize different growth phases.

  The 3 notebooks are exactly for decomposing different layers iterations. 
    - For layer (3), which is the most inner layer, [1_one_growth_curve_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/1_one_growth_curve_analysis.ipynb) is written to show how 1 single growth curve is being characterized.
    - And when we go up one more layer, for (2) and (3), [2_one_file_multiple_growth_curves_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/2_one_file_multiple_growth_curves_analysis.ipynb) is for parsing each individual growth curve from one data file.
    - And for the full feature, layer (1), (2) and (3), [3_multiple_files_analysis.ipynb](https://github.com/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/thrash-lab/sparse-growth-curve/blob/main/3_multiple_files_analysis.ipynb) reads multiple files and characterize all the growth curves within.
