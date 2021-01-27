```

  _____                      ___      _    
 |__  /___   ___  _ __ ___  / _ \    / \   
   / // _ \ / _ \| '_ ` _ \| | | |  / _ \  
  / /| (_) | (_) | | | | | | |_| | / ___ \ 
 /____\___/ \___/|_| |_| |_|\__\_\/_/   \_\
                                           

```
## Abstract
The Estimation of Model Accuracy problem is a cornerstone problem in the field of bioinformatics. Due to the resources required for X-ray crystallography and Nuclear Magnetic Resonance, computational methods for predicting a protein's tertiary structure are becoming more and more common. However, when predictions are made for proteins of which we do not know the crystal structure, we run into an issue; How do we tell how good a tertiary structure prediction is? This is the goal of the estimation of model accuracy problem. Here we introduce ZoomQA, a novel, single-model method for assessing the accuracy of a tertiary structure prediction at a residue level. ZoomQA differs from other quality assessment tools, even those that consider the 3D structure of a protein by taking this three-dimensional data one step further. ZoomQA considers the change in chemical and physical features of a fragment structure (a portion of a protein within a radius r of the target amino acid) as the radius of contact increases. ZoomQA uses fourteen physical and chemical properties of amino acids to build a comprehensive representation of every residue within a protein and grades their placement within the protein as a whole.


## Setup

#### Note, this software only works on linux environments for the time being

1. Create python virtual environment
	1. Virtualenv
		1. `pip install virtualenv` *`pip3` if you still have python2* 
		1. `python3 -m venv virtual-env-name` This creates a new virtual environment 
		1. `source virtual-env-name/bin/activate` This activates your new virtual environment 
	1. Conda 
		1. Download [Anaconda](https://www.anaconda.com/products/individual) *download the linux version to your linux machine* 
		1. Install Anaconda and follow the isntallation instructions, select yes for the init question at very end
		1. `conda activate base` to get into your 'base' environment, do not install packages to 'base' 
		1. `conda create -n virtual-env-name python=3.7`
		1. `conda activate virtual-env-name` this activates your new environment, this is where you install packages 
1. `pip install -r requirements.txt`

## Execution
1. Navigate to ZoomQA folder (currently does not work out of folder, will be updated soon)
1. `python prediction.py ./QA_examples/Input/T1096 ./TEST_OUT/`
  - This command runs the prediction and places a text file in TEST_OUT/ folder
  - An example output is provided in `Example_Out/`

#### Notes for execution
- Currently, the input data must be in a folder even if you are only running one pdb. Please put pdbs in a folder named as the target name. 
```
QA_examples
└───Input
      └───target_name
          │   input_file_1.pdb
          │   input_file_1.pdb
          │   ...
```
- Currently only works on one `target_name` as shown above, will be updated soon
