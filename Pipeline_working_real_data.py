import subprocess

Reference = input("Input reference genome: ")
subprocess.run(["samtools", "faidx", Reference], check = True)
illumina_read_1 = input("Input Illumina read 1: ")
illumina_read_2 = input("Input Illumina read 2: ")


fastq_files = [illumina_read_1, illumina_read_2]


"""
def Fastq_converter():
    for fasta_file in reads:
        fastq_file = fasta_file.rsplit('.', 1)[0] + ".fastq.gz"
        with open(fasta_file, "r") as fasta, open(fastq_file, "w") as fastq:
            read_counter=1
            for i, line in enumerate(fasta, 1):           
                seq = line.strip()
                if not seq or seq.startswith(">"):
                    continue
                header = f"read{read_counter}"
                qual = "I" * len(seq)
                fastq.write(f"@{header}\n{seq}\n+\n{qual}\n")
                read_counter+=1
        fastq_files.append(fastq_file)

Fastq_converter()    
"""    
pbam = "pipeline_mapped_reads.bam"
with open("pipeline_errors.log", "a") as errlog:
   
    p1 = subprocess.Popen(["minimap2", "-a", "-x", "sr",  Reference] + fastq_files, stdout=subprocess.PIPE, stderr=errlog)
    p2 = subprocess.Popen(["samtools", "view", "-h", "-F", "0x900", "-"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr= errlog) 
    with open(pbam, "wb") as bam_file:    
        p3 = subprocess.Popen(["samtools", "sort", "-O", "bam"], stdin=p2.stdout, stdout=bam_file, stderr=errlog)
    p3.wait()
    p1.wait()
    p2.wait()

    subprocess.run(["samtools", "index", pbam], check=True, stderr=errlog)

    with open("variant_caller_bcf.vcf", "w") as variants:
        p4 = subprocess.Popen(["bcftools", "mpileup", "-Ou", "-f", Reference, pbam], stdout=subprocess.PIPE, stderr=errlog)
        p5 = subprocess.Popen(["bcftools", "call", "-vc", "-Ov"], stdin=p4.stdout, stdout=variants, stderr=errlog)
        p4.stdout.close()
        p5.wait()
        p4.wait()
     
    subprocess.run(["bcftools", "sort", "-o", "variant_caller_bcf_sorted.vcf", "variant_caller_bcf.vcf"], check=True, stderr=errlog)
with open("pipeline_errors.log", "a") as errlog:
    
    #running snippy
    p6 = subprocess.run(["snippy", "--ref", Reference, "--outdir", "snippy_out", "--cpus", "4",  "--R1", fastq_files[0], 
                         "--R2", fastq_files[1]], check=True, stderr=errlog)

with open("pipeline_errors.log", "a") as errlog:    
    subprocess.run(["bcftools", "sort", "-O", "v", "-o", "snippy_out/snps_sorted.vcf", "snippy_out/snps.vcf"], check=True, stderr=errlog)
#setting up to use bcftools merge
with open("pipeline_errors.log", "a") as errlog:
    subprocess.run(["bgzip", "-f", "-o", "bcf_sorted.vcf.gz", "variant_caller_bcf_sorted.vcf"],
               check=True, stderr=errlog)
with open("pipeline_errors.log", "a") as errlog:
    subprocess.run(["bgzip", "-f", "-o", "snippy_sorted.vcf.gz", "snippy_out/snps_sorted.vcf"],
               check=True, stderr=errlog)


with open("new_names.txt", "w") as f:
        f.write("bcf\n")
        
with open("pipeline_errors.log", "a") as errlog:    
        
    subprocess.run(["bcftools", "reheader", "-s", "new_names.txt", "-o", "bcf_reheader.vcf.gz", "bcf_sorted.vcf.gz"],
               check=True, stderr=errlog)

with open("pipeline_errors.log", "a") as errlog:     
    subprocess.run(["bcftools", "index", "snippy_sorted.vcf.gz"], check=True, stderr=errlog)
with open("pipeline_errors.log", "a") as errlog: 
    subprocess.run(["bcftools", "index", "bcf_reheader.vcf.gz"], check=True, stderr=errlog)

with open("pipeline_errors.log", "a") as errlog:    
        
    subprocess.run(["bcftools", "merge", "--merge", "none", "snippy_sorted.vcf.gz", "bcf_reheader.vcf.gz", "-O", "z", "-o",
                    "merged.vcf.gz"],
                    check=True, stderr=errlog)
    subprocess.run(["bcftools", "index", "merged.vcf.gz"], check=True, stderr=errlog)




