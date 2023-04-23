# 0105_CBSB3_ICA
This repository is a sub-project of the CBSB3 course ICA, offered by ZJE-UOE institute. The aim is to reproduce and optimize the mathematical model of SELEX technique developed by Juli et al, 2012 via python tools (pymc3). 

**If you want to come across the source paper, please click the superlink:['Juli et al, 2012'](https://arxiv.org/abs/1205.1819)**
>**workflow for replication**
>
> 1. Extract the data from your SELEX experiment and preprocess it, including aligning the oligonucleotides and applying any necessary transformation.
> 2. Define the Gibbs free energy equation and the conditional probability distribution based on the SELEX data.
> 3. Use Monte Carlo simulation to estimate the denominator in conditional probability distribution.
> 4. Use the downhill simplex method to optimize the likelihood function of the model, which is the negative of the log-likelihood.
> 5. Extract the Gibbs free energy matrix, which represents the affinity of DNA binding sequences to the transcription factor.
> 6. Test the SELEX model with simulated and real SELEX data to validate its performance.

