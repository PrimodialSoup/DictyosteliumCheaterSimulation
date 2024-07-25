# DictyosteliumCheaterSimulation
Python model simulating the dynamics of gene-for-gene cheating and resistance in Dictyostelium Discoideum


## Usage
In a fresh conda environment (Recommended to avoid incompatibility issues): 

1. Install python
2. Install dependencies by running 'pip install -r requirements.txt'
3. Run the simulation in GUI mode with 'python DictySimulator.py'

```
conda create -n dictysim
conda activate dictysim
conda install python
pip install -r requirements.txt
python DictySimulator.py
```

## No GUI mode
To run the simulation without GUI and output directly to a .json file without visualisation, specifiy an input .json file using '--param'. For example, 'python Dictysimulator.py --param example_input.json'. An example of an acceptable input file is provided in 'example_input.json'.

## Graphing output .json files
You can visualise an output .json file using grapher.py. 'python grapher.py' will bring up a file dialog allowing you to select a .json file for visualisation.

