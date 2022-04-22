#!/bin/bash

# We assume running this from the script directory
job_directory=$PWD/job
mkdir -p ${job_directory}
src_dir="/nethome/tpan7/src/iecc/R"
data_dir="/nethome/tpan7/scgc/data/"  # need the trailing slash
out_dir_base="/nethome/tpan7/scgc/output"
log_dir="/nethome/tpan7/scgc/log"
mkdir -p $log_dir
mkdir -p $out_dir_base
cd $src_dir/..

integ="seurat"
clust="1"
f="filtered"
qc="scran"
vis="umap"
m="m"

for te in "fastwmw" "fast_t"; do # "bioqc"; do # "wilcox"; do

for prefix in "pbmc3k" "pbmc6k" "pbmc8k"; do  # 
    d=${prefix}_${f}

    for t in 1; do # 4 16 64; do # 2 8 32; do
         job_file="${job_directory}/${prefix}_${t}.job"

         nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
         mkdir -p ${out_dir_base}/${nm}

         if [[ ! -e ${log_dir}/${nm}.log ]]; then
             echo "#!/bin/bash
#SBATCH --job-name=${prefix}_${t}
#SBATCH --output=${log_dir}/log_%J.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${t}
#SBATCH --time=02:00:00
#SBATCH --partition=exclusive
Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
" > $job_file
             sbatch $job_file
         fi

## don't do no-marker
#        job_file="${job_directory}/${prefix}_${t}.job"
#
#        nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
#        mkdir -p ${out_dir_base}/${nm}
#
#        if [[ ! -e ${log_dir}/${nm}.log ]]; then
#            echo "#!/bin/bash
##SBATCH --job-name=${prefix}_${t}
##SBATCH --output=${log_dir}/log_%J.out
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=${t}
##SBATCH --time=00:10:00
##SBATCH --partition=exclusive
#Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -s &> ${log_dir}/${nm}.log
#" > $job_file
#            sbatch $job_file
#        fi

## dont do qc with marker
#         nm=${d}.${integ}.${clust}.${qc}.${te}.${vis}.${m}.${t}
#         mkdir -p ${out_dir_base}/${nm}
#         echo "#!/bin/bash
# #SBATCH --job-name=${prefix}_${t}
# #SBATCH --output=${log_dir}/log_%J.out
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=${t}
# #SBATCH --time=00:30:00
# #SBATCH --partition=exclusive
# Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -q ${qc} -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
# " > $job_file
#         sbatch $job_file

## dont do qc with no-marker
#         nm=${d}.${integ}.${clust}.${qc}.${te}.${vis}.n.${t}
#         mkdir -p ${out_dir_base}/${nm}
#         echo "#!/bin/bash
# #SBATCH --job-name=${prefix}_${t}
# #SBATCH --output=${log_dir}/log_%J.out
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=${t}
# #SBATCH --time=00:30:00
# #SBATCH --partition=exclusive
# Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -q ${qc} -t ${te} -v ${vis} -w $t -s &> ${log_dir}/${nm}.log
# " > $job_file
#         sbatch $job_file

    done
done


for prefix in "pbmc10k"; do
    d=${prefix}_${f}

    for t in 1; do # 4 16 64; do # 2 8 32; do
        job_file="${job_directory}/${prefix}_${t}.job"

        nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
        mkdir -p ${out_dir_base}/${nm}

        if [[ ! -e ${log_dir}/${nm}.log ]]; then
            echo "#!/bin/bash
#SBATCH --job-name=${prefix}_${t}
#SBATCH --output=${log_dir}/log_%J.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${t}
#SBATCH --time=04:00:00
#SBATCH --partition=exclusive
Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
" > $job_file
            sbatch $job_file
        fi
    done
done

for prefix in "pbmc33k"; do
    d=${prefix}_${f}

    for t in 1; do # 4 16 64; do #  2 8 32; do
        job_file="${job_directory}/${prefix}_${t}.job"

        nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
        mkdir -p ${out_dir_base}/${nm}

        if [[ ! -e ${log_dir}/${nm}.log ]]; then
            echo "#!/bin/bash
#SBATCH --job-name=${prefix}_${t}
#SBATCH --output=${log_dir}/log_%J.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${t}
#SBATCH --time=04:00:00
#SBATCH --partition=exclusive
Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
" > $job_file
            sbatch $job_file
        fi

    done
done


for prefix in "pbmc68k"; do
    d=${prefix}_${f}

    for t in 1; do # 4 16 64; do #  2 8 32; do
        job_file="${job_directory}/${prefix}_${t}.job"

        nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
        mkdir -p ${out_dir_base}/${nm}

        if [[ ! -e ${log_dir}/${nm}.log ]]; then
            echo "#!/bin/bash
#SBATCH --job-name=${prefix}_${t}
#SBATCH --output=${log_dir}/log_%J.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${t}
#SBATCH --time=24:00:00
#SBATCH --partition=exclusive
Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
" > $job_file
            sbatch $job_file
        fi
    done
done


for prefix in "pbmc600k"; do
    d=${prefix}_${f}

    for t in 1; do # 4 16 64; do #  2 8 32; do
        job_file="${job_directory}/${prefix}_${t}.job"

        nm=${d}.${integ}.${clust}.${f}.${te}.${vis}.${m}.${t}
        mkdir -p ${out_dir_base}/${nm}

        if [[ ! -e ${log_dir}/${nm}.log ]]; then
            echo "#!/bin/bash
#SBATCH --job-name=${prefix}_${t}
#SBATCH --output=${log_dir}/log_%J.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${t}
#SBATCH --time=24:00:00
#SBATCH --partition=exclusive
Rscript --vanilla ${src_dir}/../miscl/scgc/run_seurat.R ${data_dir} ${data_dir}${d}.csv ${out_dir_base}/${nm} ${integ} -c $clust -t ${te} -v ${vis} -w $t -m -s &> ${log_dir}/${nm}.log
" > $job_file
            sbatch $job_file
        fi
    done
done

done
