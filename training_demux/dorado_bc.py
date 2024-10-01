import os



def dorado_basecall(dorado_path, dorado_model, min_qscore, pod_dir, bam_directory, basecall_pod):
    """
    dorado_basecall takes in various inputs needed to run the Dorado basecaller
    and returns the file path to the generated BAM file 
    
    Parameters: 
    dorado_path - file pathway to Dorado
    dorado_model - canonical model to use for basecalling, refer to xr_params for more info
    xfasta_path - file pathway to an xFASTA file 
    min_qscore - minimum read quality score to allow in BAM file, set in xr_params
    pod_dir - inputted POD5 file or directory to basecall 
    bam_directory - desired output directory for outputted BAM file. 
    basecall_pod - parameter allowing for basecalling to be repeated

    Returns:
    output_bam - BAM file output file pathway
    """
    if not os.path.exists(bam_directory):
        os.makedirs(bam_directory)
    output_bam = os.path.join(bam_directory, 'bc.bam')
    
    if basecall_pod or not os.path.exists(output_bam):
        print('Xemora [STATUS] Performing basecalling using Dorado')
        #Base arguments
        dorado_args = f'--no-trim --emit-moves  {pod_dir}'

        #Optional arguments
        if min_qscore: 
            dorado_args +=f' --min-qscore {min_qscore}'

        #Dorado command 
        cmd = f'{dorado_path} basecaller {dorado_model} {dorado_args} > {output_bam}'
        os.system(cmd) 
        return output_bam
        
    else:
        print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')
        return output_bam




basecall_pod = True
dorado_path = '~/dorado-0.7.2-linux-x64/bin/dorado'
dorado_model = '~/dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0'
min_qscore = 7
pod_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/withheld_pod5'
bam_directory = '/home/xenolab/DataAnalysis/Kaplan/basecall/240930_NTC_Phusion_Withheld_Test_Set'

dorado_basecall(dorado_path, dorado_model, min_qscore, pod_dir, bam_directory, basecall_pod)
