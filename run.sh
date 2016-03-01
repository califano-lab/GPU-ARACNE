## Usage: ./bin/tester <tffile> <datafile> <ntfs> <ngenes> <nSamples> <nBootstraps> <pValue>
#make tester
expfile=$PWD/data/input_expmat_200_clean.txt
tffile=$PWD/data/tflist.txt
ngenes=`awk 'END{print NR}' $expfile`
ntfs=`awk 'END{print NR}' $tffile`
nSamples=` awk -F, 'NR==1{print NF-1}' $expfile `
nBootstraps=1
pValue=0.0000001
echo " Usage: ./bin/tester <tffile> <datafile> <ntfs> <ngenes> <nSamples> <nBootstraps> <pValue>"
cmd="$PWD/bin/runner  $tffile $expfile $ntfs $ngenes $nSamples $nBootstraps $pValue"
echo $cmd
$cmd > $PWD/data/Result_test_200.txt
