Maximum likelihood fit:
1. setup environment
   under the directory "plugins/Scripts/InitAnalysis.sh"  l2 & l11
2. change the directory of reading efficiency
   under the directory "src/Utils.cc/ReadRTEffPDF & ReadWTEffPDF"
3. choose type you want to fit
   under the directory "/python/ParameterFile.txt" 
4. compile the fitter program: make ExtractYield
5. run the fit
   eg. ./ExtractYield fit-type sample-directory eff-y/n q^2bin Maximum-order-for-eff scan-number


Method of Moment:
1. change the directory of reading efficiency
   plugins/Moment.cc l371
2. change the directory of sample
   plugins/Moment.cc l204(gen) l384(reco)
3. compile the program: make Moment
4. run 
   if gen-level: ./Moment Paramater-type gen
   if reco-level: ./Moment Paramater-type reco event-tag(mis or good)



