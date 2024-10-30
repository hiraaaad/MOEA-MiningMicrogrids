# MOEA-MiningMicrogrids
# Description

This repository contain code on "Optimal Planning of Renewable-Based Mining Microgrids: A Comparative Study of Multi-Objective Evolutionary Algorithms". The main directory include the code, results and the manuscript of the following paper.

## Usage

Main_Run_All_Algs_HPC.py". It takes 5 input variables as follows.

--nPop: population size for MOEA
--nFeval: Function evaluation
--seed_no: Seed number, int number between 0 and 30, pre-assigned for our experiments.
--DERFlag: DER combination flag, 1 for no battery and 2 for with battery cases, and 
--AlgFlag: Algorithm flag, int number between 1 and 6, these algorithms are from Pymoo.
 
Example:
```
python Main_Run_All_Algs_HPC.py --AlgFlag=1 --seed_no=0 --nPop=100 --nFeval=10000 --DERFlag=2
```
## Requirements

## Acknowledgement

## License

[MIT](https://choosealicense.com/licenses/mit/)
