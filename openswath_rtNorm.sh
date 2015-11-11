##### Step 0
# Run ExtractChromatogram for the RT Peptides
for f in ${file_path}/${file_basename}*.mzML.gz; do
  outname=`basename $f .mzML.gz`._rtnorm.chrom.mzML
  bsub -n $OMP_NUM_THREADS -R "rusage[mem=${memusage}]" ${openms_dir}OpenSwathChromatogramExtractor -is_swath -in $f -tr $irt_library \
    -out $outname -ini $ini -threads $OMP_NUM_THREADS -rt_extraction_window -1 2>&1 > /dev/null;
done

###################################
# Wait until we are done. We know that we are
# done with splitting the files because we expect the same amount of lsf files
# as we had input files.  We wait at most 3000 minutes.
for i in {0..3000}; do
    files_lsf_wc=`ls lsf* 2>/dev/null | wc `
    files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
    if [ "$files_lsf" == "$nr_files" ]; then
        echo "All files are iRT extracted"
        break;
    else 
      echo "Have ${files_lsf} files (out of ${nr_files}) so far (elapsed time ${i} minutes)"
      sleep $sleeptime;
    fi
done

files_lsf_wc=`grep 'Successfully completed' lsf* 2>/dev/null | wc `
files_lsf=`echo $files_lsf_wc | cut -f 1 -d ' ' `
if [ "$files_lsf" == "$nr_files" ]; then
    echo "All files successfully completed the iRT extraction"
else
    echo "Not all files successfully completed the iRT extraction"
    exit
fi
uptime

# merge files and run RT Normalizer
${openms_dir}FileMerger -in ${file_basename}*._rtnorm.chrom.mzML -out ${file_basename}.rtnorm.chrom.mzML
rm ${file_basename}*._rtnorm.chrom.mzML 
${openms_dir}OpenSwathRTNormalizer -in ${file_basename}.rtnorm.chrom.mzML -tr $irt_library -out ${file_basename}.rtnorm.trafoXML
tar czf runlogs_rt_normalizer.tar.gz lsf.o*
rm lsf.o*
rt_norm=${file_basename}.rtnorm.trafoXML
