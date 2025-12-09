import random

with open("Ecoli_mutated.fasta", "r") as fg:
    mutated_genome=fg.read().replace("\n", "")

print(len(mutated_genome))
#4718953

#to get 30x depth i need to multiply the length of the genome by 30 and divide by 100
print((4718953*30)/100)
#1415685.9

reads=[]
i=0
while i<1415686:
    start_pos = random.randint(0, len(mutated_genome)-100)
    read = mutated_genome[start_pos:start_pos + 100]
    reads.append((i+1, read))
    i+=1
    
print(reads[:25])
        
        

with open("Ecoli_reads.fasta", "w") as Fasta_file:
    for j, k in reads:
        Fasta_file.write(f">read_{j}\n")
        Fasta_file.write(f"{k}\n")

    

    
    
    
