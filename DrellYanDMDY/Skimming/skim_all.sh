# This file is an example of skimming and trimming.
# Before doing it, one has to carefully check all the input configs.
# The files below ARE NOT a full set

root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.data1\"\) >& log-skim-data1.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.data2\"\) >& log-skim-data2.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.data3\"\) >& log-skim-data3.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.data4\"\) >& log-skim-data4.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.data5\"\) >& log-skim-data5.txt 

root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.zttm20\"\) >& log-skim-zttm20.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.zeem20\"\) >& log-skim-zeem20.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.zeem1020\"\) >& log-skim-zeem1020.txt 
root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.zz\"\) >& log-skim-zz.txt 

root -l -b -q TrimNtuples.C+\(\"../config_files/skim_configs/trim.input.zeem20\"\) >& log-trim-zeem20.txt 
root -l -b -q TrimNtuples.C+\(\"../config_files/skim_configs/trim.input.zeem500\"\) >& log-trim-zeem500.txt 

root -l -b -q SkimNtuples.C+\(\"../config_files/skim_configs/skim.input.zeem20to500\"\) >& log-skim-zeem20to500.txt 


