#!/bin/bash

[ $# -eq 0 ] && echo "usage: $0 file (with 2cloumns samplename newname)" && exit 1

cat $1 | while read line
do
	samplename=$(echo $line|awk '{print $1}')
	newname=$(echo $line|awk '{print $2}')

	bsub -q pub.8h -n 8 -R "rusage[getdata=2,scratch=60000]" <<EOF

set -e
trap "echo JOB $samplename FAILED" EXIT

echo JOB $samplename STARTED
getdataset -v \$(searchms -m $samplename) -o \$TMPDIR -r \$TMPDIR/getdataset.out
mzXMLgz=\$(grep .mzXML.gz  \$TMPDIR/getdataset.out | awk '{print \$2}')
mzXML_to_swathmzMLgz.sh -i \$mzXMLgz -o \$TMPDIR -t 8

if [ "$newname" != "" ]
then
	mv \$TMPDIR/split_*_ms1scan.mzML.gz \$TMPDIR/${newname}_ms1.mzML.gz
	j=1
	for i in \$(ls \$TMPDIR/split_*_??.mzML.gz)
	do
  		mv -v \$i \$TMPDIR/${newname}_\$j.mzML.gz
  		j=\$(expr \$j + 1) 
	done
	[ ! -d "$newname" ] && mkdir "$newname"
	cp -v \$TMPDIR/*.mzML.gz ./$newname
else
	cp -v \$TMPDIR/*.mzML.gz ./
fi

trap - EXIT
echo JOB $samplename SUCESSFULLY FINISHED
EOF

done
