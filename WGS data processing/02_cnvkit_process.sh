if [ ${prefix} == "Normal" ];then
   normal_target_cnn=$(echo ${group_idx}|tr "," "\n"|xargs  -n 1 -I {} echo " ../Sample_{}/{}.targetcoverage.cnn ")
   ${python} ${cnvkit} gender $normal_target_cnn -o gender.txt
   ########Male
   normal_target_cnn_male=`sed '1d' gender.txt |awk '$2=="Male"'|cut -f1|tr "\n" " "`
   if echo $normal_target_cnn_male | grep -q cnn ;then
   normal_antitarget_cnn_male=${normal_target_cnn_male//targetcoverage/antitargetcoverage}
   ${python} ${cnvkit} reference $normal_target_cnn_male $normal_antitarget_cnn_male -f  ${genome} -o my_reference_male.cnn
   fi
   ########Female
   normal_target_cnn_female=`sed '1d' gender.txt |awk '$2=="Female"'|cut -f1|tr "\n" " "`
   if echo $normal_target_cnn_female | grep -q cnn; then
   normal_antitarget_cnn_female=${normal_target_cnn_female//targetcoverage/antitargetcoverage}
   ${python} ${cnvkit} reference $normal_target_cnn_female $normal_antitarget_cnn_female -f  ${genome} -o my_reference_female.cnn
   fi
fi
echo "backend:pdf" > matplotlibrc
PATH=${PATH}
gender=`grep ${ctrl_idx} ../Group_Normal/gender.txt |cut -f2`
if [ $gender = "Male" ];then
  ${python} ${cnvkit} fix ../Sample_${case_idx}/${case_idx}.targetcoverage.cnn ../Sample_${case_idx}/${case_idx}.antitargetcoverage.cnn ../Group_Normal/my_reference_male.cnn -o ${case_idx}.cnr --no-edge
  ${python} ${cnvkit} segment ${case_idx}.cnr  -t 1e-6 -o  ${case_idx}.cns -p 8
  ${python} ${cnvkit} call ${case_idx}.cns -y -v ${case_idx}_vs_${ctrl_idx}_TN.vcf -o ${case_idx}.call.cns
  ${python} ${cnvkit} gainloss ${case_idx}.cnr -y -t 0.4 -m 5 -g m -s ${case_idx}.call.cns >${case_idx}_gainloss.tsv
  ${python} ${cnvkit} scatter  ${case_idx}.cnr -s ${case_idx}.cns --title ${case_idx}_vs_${ctrl_idx} -o ${case_idx}_scatter.pdf
  ${python} ${cnvkit} diagram  ${case_idx}.cnr -s ${case_idx}.cns -o ${case_idx}_diagram.pdf
  convert -verbose -density 150 -trim ${case_idx}_scatter.pdf -quality 100 -flatten -sharpen 0x1.0 ${case_idx}_scatter.png
  convert -verbose -density 150 -trim ${case_idx}_diagram.pdf -quality 100 -flatten -sharpen 0x1.0 ${case_idx}_diagram.png
else
  ${python} ${cnvkit} fix ../Sample_${case_idx}/${case_idx}.targetcoverage.cnn ../Sample_${case_idx}/${case_idx}.antitargetcoverage.cnn ../Group_Normal/my_reference_female.cnn -o ${case_idx}.cnr --no-edge
  ${python} ${cnvkit} segment ${case_idx}.cnr -t 1e-6 -o  ${case_idx}.cns -p 8
  ${python} ${cnvkit} call ${case_idx}.cns  -v ${case_idx}_vs_${ctrl_idx}_TN.vcf -o ${case_idx}.call.cns
  ${python} ${cnvkit} gainloss ${case_idx}.cnr  -t 0.4 -m 5  -g f -s ${case_idx}.call.cns >${case_idx}_gainloss.tsv
  ${python} ${cnvkit} scatter  ${case_idx}.cnr -s ${case_idx}.cns --title ${case_idx}_vs_${ctrl_idx} -o ${case_idx}_scatter.pdf
  ${python} ${cnvkit} diagram  ${case_idx}.cnr -s ${case_idx}.cns -o ${case_idx}_diagram.pdf
  convert -verbose -density 150 -trim ${case_idx}_scatter.pdf -quality 100 -flatten -sharpen 0x1.0 ${case_idx}_scatter.png
  convert -verbose -density 150 -trim ${case_idx}_diagram.pdf -quality 100 -flatten -sharpen 0x1.0 ${case_idx}_diagram.png
fi
cns=`find ../Pair* -name "*call.cns"|xargs`
${python} ${cnvkit} heatmap $cns -o CNV_heatmap.pdf
convert -verbose -density 150 -trim CNV_heatmap.pdf -quality 100 -flatten -sharpen 0x1.0 CNV_heatmap.png

echo  -e "Pair\tGain\tLoss">CNV_gainloss.stat
for tsv in `ls ../Pair_*/*_gainloss.tsv`;do
 gain=`awk '$5>0' $tsv|wc -l`
 loss=`awk '$5<0' $tsv|wc -l`
 pair=`echo $tsv|cut -f2 -d/` 
 echo -e "$pair\t$gain\t$loss" >> CNV_gainloss.stat 
done
  
