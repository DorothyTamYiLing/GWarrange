#set default parameters in getopts
extend=100
merge=3
switch="on"

while getopts g:i:e:m:s: flag
do
    case "${flag}" in
        g) gen=${OPTARG};;
        i) IScoor=${OPTARG};;
        e) extend=${OPTARG};;
        m) merge=${OPTARG};;
        s) switch=${OPTARG};;
    esac
done
echo "gen: $gen";
echo "IS blast output: $IScoor"
echo "IS extend: $extend";
echo "IS merge: $merge";
echo "switch for turning on/off minimal ISmerge: $switch";

#merge and extending IS
Rscript ./scripts/merge_IS.R --input $IScoor --extend $extend --merge $merge 

python3 ./scripts/iSreplace_2col.py --input ${gen}  --coor ext${extend}_merge${merge}_mergedIS.txt --out ext${extend}_merge${merge}_ISreplaced_genomes

if [[ ($merge -gt 3 || $extend -gt 100) && $switch == "on" ]]  #if merge or extend arguments is provided and larger than default values, and switch for minimal ISmerg is on as default
then
echo "also perform merging overlapping IS only"

Rscript ./scripts/merge_IS.R --input $IScoor --extend 100 --merge 3

python3 ./scripts/iSreplace_2col.py --input ${gen}  --coor ext100_merge3_mergedIS.txt --out ext100_merge3_ISreplaced_genomes
fi
