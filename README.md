# DataMCreweighting
Reweighting code of data/MC difference



## Apply the existing model without comparing


Adapt existing model 

 ***fulldata_trained_model_2018_gitv2_onlyRECO_addPhi.json***

to other q2Bins you can just use 

***Application.py***

Change the input MC file in line 515, input json file in line 577, and output MC weights tree in line 595

Then use 

```sh
python Application.py q2Bin parity year
```

If need to save the weights in ntuple

```sh
python Addreweight.py q2Bin parity year
```


If you want to train from beginning, please follow 

## Create Dataset

```sh
root -l 'createDataset.cc(year,q2Bin,data,deno,num)'
```
Then we can generate part of data/MC ntuple after all selections and split trk variables in Pt (trk1 indicates the leading trk, trk2 indicates subleading track) (Split several parts due to limitations of ntuple size)

Then add sWeights to data ntuple, merge all files and divide into two parities

```sh
python AddSweight.py
```
For MC we just need to merge all root files

```sh
python MCmerge.py q2Bin parity year
```

## Generate model

We just use even parity Jpsi channel to generate reweighting model

```sh
python DoReweight_readindat_fulldatatraining.py q2Bin parity year
```
Here we use "even" parity Jpsi channel events to train the model
## Apply model

Then we can apply model to each q2Bin (here use odd parity events) 

```sh
python Application_DoReweight_readindat_fulldatatraining.py q2Bin parity year
```
remember that we can only compare data and MC distributions in Jpsi channel and Psi(2S) channel

## Generate final ntuple

Finally we merge data/MC weight values with original ntuple

```sh
python Addreweight.py q2Bin parity year

