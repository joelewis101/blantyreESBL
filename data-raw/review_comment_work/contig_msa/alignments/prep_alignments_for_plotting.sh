#!/bin/zsh

# Change all .'s to N's for PopGenome R package

FILES_ARRAY=($(ls -1 *.aln ))
for CURRENTREF in "${FILES_ARRAY[@]}"
do
echo "Current file:" $CURRENTREF
OUTPREFIX=$(echo $CURRENTREF | rev | cut -c5- | rev)
echo "Making directory..." $OUTPREFIX
rm -R $OUTPREFIX
mkdir $OUTPREFIX
OUTFILE=$OUTPREFIX.mod.aln
echo "Changing all .s to Ns, writing to" $OUTPREFIX/$OUTFILE "..."
sed -E 's/\.|\*/N/g' $CURRENTREF > $OUTPREFIX/$OUTFILE
echo "Done"
done
