
#Convert to unified format
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo Creatinine | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo CRP | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo Glucose | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo HDL | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo HGB | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo HOMA_b | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo HOMA_ir | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo IL6 | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo Insulin | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo MCHC | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo MCH | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo MCV | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo PCV | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo PLT | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo RBC | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo TC | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo TG | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo UricAcid | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo WBC | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"
bsub -G team170 -n1 -R "span[hosts=1] select[mem>10000] rusage[mem=10000]" -q normal -M 10000 -o unify_gwas.%J.jobout "echo LDL | python ~/software/utils/GWAS/metabolicToUnified.py --indir Metabolic_GWAS/ --outdir Inflammatory_GWAS/unsorted/"

#Sort an index summary stats
echo Creatinine.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo CRP.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo Glucose.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo HDL.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo HGB.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo HOMA_b.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo HOMA_ir.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo IL6.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo Insulin.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo MCHC.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo MCH.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo MCV.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo PCV.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo PLT.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo RBC.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo TC.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo TG.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo UricAcid.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo WBC.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
echo LDL.unified | python ~/software/utils/GWAS/tabixGWASSummaryStats.py --indir datasets/Inflammatory_GWAS/unsorted/
