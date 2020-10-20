# ZoomQA

## Abstract
The Estimation of Model Accuracy problem is a cornerstone problem in the field of bioinformatics. Due to the resources required for X-ray crystallography and Nuclear Magnetic Resonance, computational methods for predicting a protein's tertiary structure are becoming more and more common. However, when predictions are made for proteins of which we do not know the crystal structure, we run into an issue; How do we tell how good a tertiary structure prediction is? This is the goal of the estimation of model accuracy problem. Here we introduce ZoomQA, a novel, single-model method for assessing the accuracy of a tertiary structure prediction at a residue level. ZoomQA differs from other quality assessment tools, even those that consider the 3D structure of a protein by taking this three-dimensional data one step further. ZoomQA considers the change in chemical and physical features of a fragment structure (a portion of a protein within a radius r of the target amino acid) as the radius of contact increases. ZoomQA uses fourteen physical and chemical properties of amino acids to build a comprehensive representation of every residue within a protein and grades their placement within the protein as a whole.


## Setup

#### Note, this software only works on linux environments for the time being

1. Create python virtual environment
  1. Virtualenv
    1. `pip install virtualenv` *`pip3` if you still have python2*
    1. `python3 -m venv virtual-env-name`
    1. `source virtual-env-name/bin/activate`
  1. Conda *after downloading, and activating base*
    1. `conda create -n virtual-env-name python=3.7`
    1. `conda activate virtual-env-name`
1. `pip install -r requirements.txt`

## Execution
1. Navigate to ZoomQA folder (currently does not work out of folder, will be updated soon)
1. `python prediction.py ./QA_examples/Input/ ./TEST_OUT/`
  - This command runs the prediction and places a text file in TEST_OUT/ folder
  - An example output is provided in `Example_Out/`

#### Notes for execution
- Currently, the input data must be nested in a folder, example below 
