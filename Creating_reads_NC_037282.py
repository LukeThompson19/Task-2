import random

with open("NC_037282_mutated.fasta", "r") as fg:
    mutated_genome=fg.read().replace("\n", "")

print(len(mutated_genome))
#2072258

#to get 30x depth i need to multiply the length of the genome by 30 and divide by 100
print((2072258*30)/100)
#621677.4

reads=[]
i=0
while i<621678:
    start_pos = random.randint(0, len(mutated_genome)-100)
    read = mutated_genome[start_pos:start_pos + 100]
    reads.append((i+1, read))
    i+=1
    
print(reads[:25])
        
        

with open("NC_037282_reads.fasta", "w") as Fasta_file:
    for j, k in reads:
        Fasta_file.write(f">read_{j}\n")
        Fasta_file.write(f"{k}\n")

    

    
    
    
