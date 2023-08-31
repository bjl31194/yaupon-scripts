#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=HAPevery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=168:00:00
#SBATCH --output=GATKhapevery.out
#SBATCH --error=GATK_SVevery.err

ml GATK/4.2.5.0-GCCcore-8.3.0-Java-1.8


#gatk CombineGVCFs -R ZJN_r1.1.fa --variant 13-10Hap.g.vcf	--variant 28-11Hap.g.vcf	--variant 36-1Hap.g.vcf	--variant 46-6Hap.g.vcf	--variant 56-3LPHap.g.vcf	--variant 67-8LPHap.g.vcf	--variant 82-12Hap.g.vcf	--variant 92-8LPHap.g.vcf	--variant t4-1Hap.g.vcf	--variant 13-15Hap.g.vcf	--variant 28-12Hap.g.vcf	--variant 36-3Hap.g.vcf	--variant 46-7Hap.g.vcf	--variant 56-8LPHap.g.vcf	--variant 68-10Hap.g.vcf	--variant 82-12LPHap.g.vcf	--variant 93-13Hap.g.vcf	--variant T43-15Hap.g.vcf	--variant 13-7LPHap.g.vcf	--variant 28-6Hap.g.vcf	--variant 36-4Hap.g.vcf	--variant 46-8Hap.g.vcf	--variant 57-10LPHap.g.vcf	--variant 68-2LPHap.g.vcf	--variant 82-15Hap.g.vcf	--variant 93-13LPHap.g.vcf	--variant T44-12Hap.g.vcf	--variant 13-8Hap.g.vcf	--variant 28-8Hap.g.vcf	--variant 36-5LPHap.g.vcf	--variant 46-8LPHap.g.vcf	--variant 57-13LPHap.g.vcf	--variant 68-4Hap.g.vcf	--variant 82-16Hap.g.vcf	--variant 93-14LPHap.g.vcf	--variant t4-7Hap.g.vcf	--variant 13-9Hap.g.vcf	--variant 29-10Hap.g.vcf	--variant 36-7Hap.g.vcf	--variant 46-9Hap.g.vcf	--variant 57-16LPHap.g.vcf	--variant 68-6LPHap.g.vcf	--variant 82-4LPHap.g.vcf	--variant 93-16Hap.g.vcf	--variant T49-7Hap.g.vcf	--variant 14-12Hap.g.vcf	--variant 29-11LPHap.g.vcf	--variant 36-9Hap.g.vcf	--variant 47-11Hap.g.vcf	--variant 57-6LPHap.g.vcf	--variant 68-7LPHap.g.vcf	--variant 82-6Hap.g.vcf	--variant 94-14LPHap.g.vcf	--variant t5-11Hap.g.vcf	--variant 15-10Hap.g.vcf	--variant 29-12LPHap.g.vcf	--variant 37-13LPHap.g.vcf	--variant 47-12Hap.g.vcf	--variant 58-13Hap.g.vcf	--variant 69-14Hap.g.vcf	--variant 82-8Hap.g.vcf	--variant 94-16LPHap.g.vcf	--variant T51-3Hap.g.vcf	--variant 15-13Hap.g.vcf	--variant 29-14LPHap.g.vcf	--variant 37-15Hap.g.vcf	--variant 47-16Hap.g.vcf	--variant 58-14Hap.g.vcf	--variant 69-15Hap.g.vcf	--variant 83-11Hap.g.vcf	--variant 94-3Hap.g.vcf	--variant t5-16Hap.g.vcf	--variant 15-3Hap.g.vcf	--variant 29-15Hap.g.vcf	--variant 37-2Hap.g.vcf	--variant 47-1Hap.g.vcf	--variant 58-15Hap.g.vcf	--variant 70-16LPHap.g.vcf	--variant 83-11LPHap.g.vcf	--variant 94-6Hap.g.vcf	--variant T5-16PPHap.g.vcf	--variant 15-4Hap.g.vcf	--variant 29-16Hap.g.vcf	--variant 37-3Hap.g.vcf	--variant 47-2Hap.g.vcf	--variant 58-1LPHap.g.vcf	--variant 70-1LPHap.g.vcf	--variant 83-13Hap.g.vcf	--variant 94-7Hap.g.vcf	--variant t5-2Hap.g.vcf	--variant 15-5Hap.g.vcf	--variant 29-2Hap.g.vcf	--variant 37-6Hap.g.vcf	--variant 47-3Hap.g.vcf	--variant 58-2LPHap.g.vcf	--variant 70-4LPHap.g.vcf	--variant 83-15LPHap.g.vcf	--variant 95-11Hap.g.vcf	--variant t5-7Hap.g.vcf	--variant 15-6Hap.g.vcf	--variant 29-4Hap.g.vcf	--variant 37-7Hap.g.vcf	--variant 47-4Hap.g.vcf	--variant 58-5LPHap.g.vcf	--variant 70-5Hap.g.vcf	--variant 83-16Hap.g.vcf	--variant 95-13Hap.g.vcf	--variant T58-3Hap.g.vcf	--variant 15-9Hap.g.vcf	--variant 29-5Hap.g.vcf	--variant 37-8Hap.g.vcf	--variant 47-5Hap.g.vcf	--variant 58-6Hap.g.vcf	--variant 70-6LPHap.g.vcf	--variant 83-6Hap.g.vcf	--variant 95-16Hap.g.vcf	--variant T60-2Hap.g.vcf	--variant 16-4LPHap.g.vcf	--variant 29-6Hap.g.vcf	--variant 38-12Hap.g.vcf	--variant 47-7Hap.g.vcf	--variant 59-11Hap.g.vcf	--variant 70-8LPHap.g.vcf	--variant 84-5Hap.g.vcf	--variant 95-3Hap.g.vcf	--variant T6-1Hap.g.vcf	--variant 16-8Hap.g.vcf	--variant 29-8LPHap.g.vcf	--variant 38-15Hap.g.vcf	--variant 47-8Hap.g.vcf	--variant 59-12Hap.g.vcf	--variant 70-9LPHap.g.vcf	--variant 84-7Hap.g.vcf	--variant 95-4Hap.g.vcf	--variant T62-9Hap.g.vcf	--variant 17-11Hap.g.vcf	--variant 30-10LPHap.g.vcf	--variant 38-16Hap.g.vcf	--variant 47-9Hap.g.vcf	--variant 59-13Hap.g.vcf	--variant 71-12Hap.g.vcf	--variant 84-7LPHap.g.vcf	--variant 95-6Hap.g.vcf	--variant t6-29LPHap.g.vcf	--variant 17-15Hap.g.vcf	--variant 30-12LPHap.g.vcf	--variant 38-1Hap.g.vcf	--variant 48-13Hap.g.vcf	--variant 59-2LPHap.g.vcf	--variant 71-2LPHap.g.vcf	--variant 84-8Hap.g.vcf	--variant F1FirstHap.g.vcf	--variant t6-3Hap.g.vcf	--variant 17-16LPHap.g.vcf	--variant 30-13LPHap.g.vcf	--variant 38-4Hap.g.vcf	--variant 48-2Hap.g.vcf	--variant 59-6LPHap.g.vcf	--variant 71-5LPHap.g.vcf	--variant 85-14Hap.g.vcf	--variant F1LPHap.g.vcf	--variant t6-4LPHap.g.vcf	--variant 17-3Hap.g.vcf	--variant 30-15LPHap.g.vcf	--variant 38-5Hap.g.vcf	--variant 48-3Hap.g.vcf	--variant 59-8LPHap.g.vcf	--variant 71-9LPHap.g.vcf	--variant 85-16Hap.g.vcf	--variant MeyerFirstHap.g.vcf	--variant T65-13Hap.g.vcf	--variant 17-6Hap.g.vcf	--variant 30-16LPHap.g.vcf	--variant 38-7Hap.g.vcf	--variant 48-3LPHap.g.vcf	--variant 59-9LPHap.g.vcf	--variant 72-11Hap.g.vcf	--variant 85-16LPHap.g.vcf	--variant MeyerLPHap.g.vcf	--variant T66-9Hap.g.vcf	--variant 17-8Hap.g.vcf	--variant 30-1Hap.g.vcf	--variant 38-8LPHap.g.vcf	--variant 48-5Hap.g.vcf	--variant 60-11LPHap.g.vcf	--variant 72-13Hap.g.vcf	--variant 85-3Hap.g.vcf	--variant Mraw_SNPs1.vcf	--variant T67-5Hap.g.vcf	--variant 18-15Hap.g.vcf	--variant 30-6Hap.g.vcf	--variant 39-10Hap.g.vcf	--variant 49-11Hap.g.vcf	--variant 60-14LPHap.g.vcf	--variant 72-16LPHap.g.vcf	--variant 85-6Hap.g.vcf	--variant Mraw_SNPs.vcf	--variant t6-7LPHap.g.vcf	--variant 18-16Hap.g.vcf	--variant 30-6LPHap.g.vcf	--variant 39-12Hap.g.vcf	--variant 49-12Hap.g.vcf	--variant 60-16Hap.g.vcf	--variant 72-8Hap.g.vcf	--variant 85-6LPHap.g.vcf	--variant PI231146First2Hap.g.vcf	--variant T68-7.g.gatkv3.vcf	--variant 18-2Hap.g.vcf	--variant 30-7Hap.g.vcf	--variant 39-13Hap.g.vcf	--variant 49-14Hap.g.vcf	--variant 60-1LPHap.g.vcf	--variant 72-9Hap.g.vcf	--variant 86-3Hap.g.vcf	--variant PI231146LPHap.g.vcf	--variant T68-7Hap.g.vcf	--variant 19-12LPHap.g.vcf	--variant 30-8Hap.g.vcf	--variant 39-14Hap.g.vcf	--variant 49-16Hap.g.vcf	--variant 60-2Hap.g.vcf	--variant 73-10LPHap.g.vcf	--variant 86-4Hap.g.vcf	--variant PIFirstHap.g.vcf	--variant T68-7test2parameters.g.vcf	--variant 19-6Hap.g.vcf	--variant 30-9Hap.g.vcf	--variant 39-4Hap.g.vcf	--variant 49-3Hap.g.vcf	--variant 60-4LPHap.g.vcf	--variant 73-12Hap.g.vcf	--variant 86-7LPHap.g.vcf	--variant t10-10Hap.g.vcf	--variant T68-7testv4NoP.g.vcf	--variant 19-9Hap.g.vcf	--variant 31-13LPHap.g.vcf	--variant 39-6Hap.g.vcf	--variant 49-7Hap.g.vcf	--variant 60-8-aLPHap.g.vcf	--variant 73-15LPHap.g.vcf	--variant 87-11Hap.g.vcf	--variant t10-10LPHap.g.vcf	--variant T68-7testv4SM.g.vcf	--variant 20-10Hap.g.vcf	--variant 31-14aLPHap.g.vcf	--variant 39-7Hap.g.vcf	--variant 50-11Hap.g.vcf	--variant 60-8LPHap.g.vcf	--variant 73-2LPHap.g.vcf	--variant 87-12Hap.g.vcf	--variant t10-3Hap.g.vcf	--variant t7-12LPHap.g.vcf	--variant 20-11Hap.g.vcf	--variant 31-14LPHap.g.vcf	--variant 40-12Hap.g.vcf	--variant 50-14Hap.g.vcf	--variant 60-9Hap.g.vcf	--variant 73-3Hap.g.vcf	--variant 87-13Hap.g.vcf	--variant t10-5Hap.g.vcf	--variant t7-16LPHap.g.vcf	--variant 20-14Hap.g.vcf	--variant 31-15Hap.g.vcf	--variant 40-13Hap.g.vcf	--variant 50-15Hap.g.vcf	--variant 61-10LPHap.g.vcf	--variant 73-4Hap.g.vcf	--variant 87-14Hap.g.vcf	--variant t11-13Hap.g.vcf	--variant T73-10Hap.g.vcf	--variant 20-2Hap.g.vcf	--variant 31-1LPHap.g.vcf	--variant 40-15Hap.g.vcf	--variant 50-1Hap.g.vcf	--variant 61-11LPHap.g.vcf	--variant 73-5Hap.g.vcf	--variant 87-14LPHap.g.vcf	--variant t11-2Hap.g.vcf	--variant t7-3Hap.g.vcf	--variant 20-6Hap.g.vcf	--variant 31-2Hap.g.vcf	--variant 40-3Hap.g.vcf	--variant 50-3Hap.g.vcf	--variant 61-15LPHap.g.vcf	--variant 74-10LPHap.g.vcf	--variant 87-15LPHap.g.vcf	--variant t11-7Hap.g.vcf	--variant T75-8Hap.g.vcf	--variant 20-8Hap.g.vcf	--variant 31-7LPHap.g.vcf	--variant 40-5Hap.g.vcf	--variant 50-5Hap.g.vcf	--variant 61-16LPHap.g.vcf	--variant 74-14LPHap.g.vcf	--variant 87-4Hap.g.vcf	--variant t11-8Hap.g.vcf	--variant t7-5Hap.g.vcf	--variant 21-12Hap.g.vcf	--variant 31-8LPHap.g.vcf	--variant 40-7Hap.g.vcf	--variant 51-12Hap.g.vcf	--variant 61-2LPHap.g.vcf	--variant 76-2Hap.g.vcf	--variant 87-4LPHap.g.vcf	--variant t11-9Hap.g.vcf	--variant t7-6Hap.g.vcf	--variant 21-13Hap.g.vcf	--variant 32-10Hap.g.vcf	--variant 40-8LPHap.g.vcf	--variant 51-16Hap.g.vcf	--variant 61-6LPHap.g.vcf	--variant 76-4Hap.g.vcf	--variant 87-7Hap.g.vcf	--variant t12-10LPHap.g.vcf	--variant T77-13Hap.g.vcf	--variant 21-14Hap.g.vcf	--variant 32-11LPHap.g.vcf	--variant 40-9Hap.g.vcf	--variant 51-1Hap.g.vcf	--variant 61-7LPHap.g.vcf	--variant 76-6Hap.g.vcf	--variant 87-9LPHap.g.vcf	--variant t12-11Hap.g.vcf	--variant t7-7Hap.g.vcf	--variant 21-7Hap.g.vcf	--variant 32-13Hap.g.vcf	--variant 41-11Hap.g.vcf	--variant 51-2Hap.g.vcf	--variant 62-11LPHap.g.vcf	--variant 77-12Hap.g.vcf	--variant 88-12Hap.g.vcf	--variant T12-14Hap.g.vcf	--variant T8-10Hap.g.vcf	--variant 22-10Hap.g.vcf	--variant 32-16Hap.g.vcf	--variant 41-12Hap.g.vcf	--variant 51-4Hap.g.vcf	--variant 62-14LPHap.g.vcf	--variant 77-12LPHap.g.vcf	--variant 88-13Hap.g.vcf	--variant t12-15Hap.g.vcf	--variant t8-10LPHap.g.vcf	--variant 22-11Hap.g.vcf	--variant 32-2LPHap.g.vcf	--variant 41-14Hap.g.vcf	--variant 51-7Hap.g.vcf	--variant 62-2LPHap.g.vcf	--variant 77-15Hap.g.vcf	--variant 88-15Hap.g.vcf	--variant t12-2Hap.g.vcf	--variant t8-12Hap.g.vcf	--variant 22-14Hap.g.vcf	--variant 32-3Hap.g.vcf	--variant 41-5Hap.g.vcf	--variant 51-8Hap.g.vcf	--variant 62-3LPHap.g.vcf	--variant 77-1Hap.g.vcf	--variant 88-2Hap.g.vcf	--variant T15-16Hap.g.vcf	--variant t8-14Hap.g.vcf	--variant 22-15Hap.g.vcf	--variant 32-8Hap.g.vcf	--variant 41-6Hap.g.vcf	--variant 51-8LPHap.g.vcf	--variant 62-4LPHap.g.vcf	--variant 77-4Hap.g.vcf	--variant 88-3Hap.g.vcf	--variant t15-1Hap.g.vcf	--variant t8-15Hap.g.vcf	--variant 22-6LPHap.g.vcf	--variant 32-9Hap.g.vcf	--variant 41-8Hap.g.vcf	--variant 52-10Hap.g.vcf	--variant 62-5LPHap.g.vcf	--variant 77-7Hap.g.vcf	--variant 88-7Hap.g.vcf	--variant t15-2Hap.g.vcf	--variant t8-15LPHap.g.vcf	--variant 22-8Hap.g.vcf	--variant 33-11Hap.g.vcf	--variant 41-9Hap.g.vcf	--variant 52-14Hap.g.vcf	--variant 62-6LPHap.g.vcf	--variant 77-7LPHap.g.vcf	--variant 89-12Hap.g.vcf	--variant t1-5Hap.g.vcf	--variant t8-19LPHap.g.vcf	--variant 22-9LPHap.g.vcf	--variant 33-14Hap.g.vcf	--variant 42-5Hap.g.vcf	--variant 52-14LPHap.g.vcf	--variant 62-7LPHap.g.vcf	--variant 77-9LPHap.g.vcf	--variant 89-13Hap.g.vcf	--variant T19-12Hap.g.vcf	--variant t8-2Hap.g.vcf	--variant 23-12Hap.g.vcf	--variant 33-15Hap.g.vcf	--variant 43-11Hap.g.vcf	--variant 52-1Hap.g.vcf	--variant 62-8LPHap.g.vcf	--variant 78-12Hap.g.vcf	--variant 89-2Hap.g.vcf	--variant t1-9Hap.g.vcf	--variant T83-11Hap.g.vcf	--variant 23-15Hap.g.vcf	--variant 33-1Hap.g.vcf	--variant 43-14LPHap.g.vcf	--variant 52-3Hap.g.vcf	--variant 63-1Hap.g.vcf	--variant 78-16Hap.g.vcf	--variant 89-4Hap.g.vcf	--variant T20-10Hap.g.vcf	--variant t8-3Hap.g.vcf	--variant 23-16LPHap.g.vcf	--variant 33-2Hap.g.vcf	--variant 43-15Hap.g.vcf	--variant 52-5Hap.g.vcf	--variant 63-1LPHap.g.vcf	--variant 78-16LPHap.g.vcf	--variant 89-9Hap.g.vcf	--variant t2-10Hap.g.vcf	--variant t8-4Hap.g.vcf	--variant 23-2LPHap.g.vcf	--variant 33-3Hap.g.vcf	--variant 43-16Hap.g.vcf	--variant 53-12Hap.g.vcf	--variant 63-3LPHap.g.vcf	--variant 78-5Hap.g.vcf	--variant 90-14LPHap.g.vcf	--variant t2-13Hap.g.vcf	--variant T87-13Hap.g.vcf	--variant 23-3Hap.g.vcf	--variant 33-4Hap.g.vcf	--variant 43-3Hap.g.vcf	--variant 53-14Hap.g.vcf	--variant 63-5LPHap.g.vcf	--variant 78-6LPHap.g.vcf	--variant 90-2Hap.g.vcf	--variant t2-1Hap.g.vcf	--variant T87-15Hap.g.vcf	--variant 23-6Hap.g.vcf	--variant 33-6LPHap.g.vcf	--variant 43-4Hap.g.vcf	--variant 53-16Hap.g.vcf	--variant 63-6LPHap.g.vcf	--variant 78-8Hap.g.vcf	--variant 90-2LPHap.g.vcf	--variant T22-9Hap.g.vcf	--variant t8-7Hap.g.vcf	--variant 23-8LPHap.g.vcf	--variant 33-9Hap.g.vcf	--variant 44-13Hap.g.vcf	--variant 53-1Hap.g.vcf	--variant 64-11Hap.g.vcf	--variant 78-9LPHap.g.vcf	--variant 90-3Hap.g.vcf	--variant T23-2Hap.g.vcf	--variant t8-7LPHap.g.vcf	--variant 23-9Hap.g.vcf	--variant 34-10Hap.g.vcf	--variant 44-15Hap.g.vcf	--variant 53-2Hap.g.vcf	--variant 64-12Hap.g.vcf	--variant 79-13Hap.g.vcf	--variant 90-4Hap.g.vcf	--variant T25-6Hap.g.vcf	--variant t8-8Hap.g.vcf	--variant 24-10LPHap.g.vcf	--variant 34-13Hap.g.vcf	--variant 44-2LPHap.g.vcf	--variant 53-4LPHap.g.vcf	--variant 64-13Hap.g.vcf	--variant 79-14Hap.g.vcf	--variant 90-5Hap.g.vcf	--variant t2-5Hap.g.vcf	--variant T90-2Hap.g.vcf	--variant 24-1Hap.g.vcf	--variant 34-16Hap.g.vcf	--variant 44-3Hap.g.vcf	--variant 53-5Hap.g.vcf	--variant 64-4LPHap.g.vcf	--variant 79-2Hap.g.vcf	--variant 90-6Hap.g.vcf	--variant T27-4Hap.g.vcf	--variant t9-12LPHap.g.vcf	--variant 24-9Hap.g.vcf	--variant 34-5Hap.g.vcf	--variant 44-3LPHap.g.vcf	--variant 53-8Hap.g.vcf	--variant 64-6Hap.g.vcf	--variant 79-3Hap.g.vcf	--variant 90-7Hap.g.vcf	--variant t2-7LPHap.g.vcf	--variant t9-13Hap.g.vcf	--variant 25-15Hap.g.vcf	--variant 34-7Hap.g.vcf	--variant 44-6Hap.g.vcf	--variant 54-11Hap.g.vcf	--variant 64-7LPHap.g.vcf	--variant 79-5Hap.g.vcf	--variant 91-11Hap.g.vcf	--variant T30-11Hap.g.vcf	--variant T91-3Hap.g.vcf	--variant 25-2Hap.g.vcf	--variant 34-8Hap.g.vcf	--variant 44-7Hap.g.vcf	--variant 54-14Hap.g.vcf	--variant 64-8LPHap.g.vcf	--variant 79-7Hap.g.vcf	--variant 91-14Hap.g.vcf	--variant t3-13LPHap.g.vcf	--variant t9-16Hap.g.vcf	--variant 26-12Hap.g.vcf	--variant 34-9Hap.g.vcf	--variant 44-8Hap.g.vcf	--variant 54-16LPHap.g.vcf	--variant 64-9LPHap.g.vcf	--variant 79-8Hap.g.vcf	--variant 91-2Hap.g.vcf	--variant t3-14LPHap.g.vcf	--variant t9-2Hap.g.vcf	--variant 26-13Hap.g.vcf	--variant 35-10LPHap.g.vcf	--variant 44-9Hap.g.vcf	--variant 54-4Hap.g.vcf	--variant 66-13LPHap.g.vcf	--variant 80-15Hap.g.vcf	--variant 91-3LPHap.g.vcf	--variant t3-15aLPHap.g.vcf	--variant t9-3Hap.g.vcf	--variant 26-15Hap.g.vcf	--variant 35-16LPHap.g.vcf	--variant 45-11Hap.g.vcf	--variant 54-6LPHap.g.vcf	--variant 66-14LPHap.g.vcf	--variant 80-16Hap.g.vcf	--variant 91-4Hap.g.vcf	--variant t3-15LPHap.g.vcf	--variant t9-3LPHap.g.vcf	--variant 26-1Hap.g.vcf	--variant 35-4Hap.g.vcf	--variant 45-14Hap.g.vcf	--variant 54-7LPHap.g.vcf	--variant 66-7LPHap.g.vcf	--variant 80-1Hap.g.vcf	--variant 91-4LPHap.g.vcf	--variant t3-16Hap.g.vcf	--variant t9-4Hap.g.vcf	--variant 35-6Hap.g.vcf	--variant 45-1Hap.g.vcf	--variant 54-8LPHap.g.vcf	--variant 66-8Hap.g.vcf	--variant 80-7Hap.g.vcf	--variant 91-7Hap.g.vcf	--variant T3-1Hap.g.vcf	--variant t9-7Hap.g.vcf	--variant 26-9LPHap.g.vcf	--variant 35-7Hap.g.vcf	--variant 45-2Hap.g.vcf	--variant 54-9LPHap.g.vcf	--variant 66-9LPHap.g.vcf	--variant 80-9Hap.g.vcf	--variant 91-8Hap.g.vcf	--variant T34-8Hap.g.vcf	--variant t9-8LPHap.g.vcf	--variant 27-10Hap.g.vcf	--variant 35-8Hap.g.vcf	--variant 45-7Hap.g.vcf	--variant 55-10LPHap.g.vcf	--variant 67-10LPHap.g.vcf	--variant 81-3Hap.g.vcf	--variant 91-9Hap.g.vcf	--variant T35-16Hap.g.vcf	--variant 27-13Hap.g.vcf	--variant 35-9Hap.g.vcf	--variant 45-8Hap.g.vcf	--variant 55-11LPHap.g.vcf	--variant 67-11LPHap.g.vcf	--variant 81-6LPHap.g.vcf	--variant 91-9LPHap.g.vcf	--variant T36-11Hap.g.vcf	--variant 27-14Hap.g.vcf	--variant 36-10Hap.g.vcf	--variant 46-10Hap.g.vcf	--variant 55-16Hap.g.vcf	--variant 67-13LPHap.g.vcf	--variant 81-7LPHap.g.vcf	--variant 92-11Hap.g.vcf	--variant t3-8Hap.g.vcf	--variant 27-16Hap.g.vcf	--variant 36-11Hap.g.vcf	--variant 46-10LPHap.g.vcf	--variant 55-1LPHap.g.vcf	--variant 67-16LPHap.g.vcf	--variant 81-9Hap.g.vcf	--variant 92-13Hap.g.vcf	--variant t3-9Hap.g.vcf	--variant 27-3LPHap.g.vcf	--variant 36-12Hap.g.vcf	--variant 46-16Hap.g.vcf	--variant 55-4LPHap.g.vcf	--variant 67-3Hap.g.vcf	--variant 81-9LPHap.g.vcf	--variant 92-15Hap.g.vcf	--variant t4-10Hap.g.vcf	--variant 27-7Hap.g.vcf	-

cd /scratch/sp27971/Zoysia/UnifiedGenotyper/HaplotypeCaller

gatk CombineGVCFs -R ZJN_r1.1.fa --variant gVCF.list -O HapCombineGVCF.g.vcf

gatk GenotypeGVCFs -R ZJN_r1.1.fa -O HapGenotypeGVCF.vcf --V HapCombineGVCF.g.vcf

gatk SelectVariants -R ZJN_r1.1.fa -O SelectVariantfromHap.vcf --select-type-to-include SNP -V HapGenotypeGVCF.vcf