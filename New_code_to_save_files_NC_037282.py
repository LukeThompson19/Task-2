import random
import copy
import json
genome = "NC_037282.1.fasta"

header = ""
sequence = []


bases = ["A", "T", "G", "C"]
#get the sequence and store the header seperately
with open(genome, "r") as file:
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            header = line
            continue
        sequence.extend(list(line))
            
#print(sequence[:10])

#snp position generator
snp_position= []

occupied_snp_position = set()
def mutater_snp_pos(sequence_to_mutate_snps):
    i=0
    while i < 300:
        snp_pos=random.randint(0, len(sequence_to_mutate_snps)-1)
        if snp_pos in occupied_snp_position:
            continue
        snp_position.append(snp_pos)
        occupied_snp_position.add(snp_pos)
        i=i+1
    
        
mutater_snp_pos(sequence)        

#print(len(snp_position))

snp_list =[]

#snp mutater
def snp_converter(sequence_to_add_snps):
    for i in snp_position:
        original_base = sequence_to_add_snps[i]
        mutant_base = random.choice([b for b in bases if b != original_base])       
        snp_list.append(("snp", i, i, [original_base], [mutant_base]))
        
    return snp_list

snp_converter(sequence)
#print(snp_list)

#add to dictionary
variant_dict = []

keys = ["mutation_type", "start", "end", "original", "new"]

for i in snp_list:
    variant_dict.append(dict(zip(keys, i)))
    







insertions = []
occupied_insertion_range=set()
#insert position generator
def mutater_insertion_pos(sequence_to_mutate_insertion):
    i=0
    while i < 10:
        
        insert_pos=random.randint(0, len(sequence_to_mutate_insertion)-1)
        
        
        if insert_pos in snp_position:
            continue
        if insert_pos in occupied_insertion_range:
            continue
        insertions.append(insert_pos)
        occupied_insertion_range.add(insert_pos)
        i=i+1
    return insertions

mutater_insertion_pos(sequence)
#print(insertions)

#insertion mutater

insert_list = []
def insertion_converter(sequence_to_add_insertions):
    for start in insertions:       
        insert_length = random.randint(1, 10)  
        new_seq = [random.choice(bases) for _ in range(insert_length)] 
        insert_list.append(("insert", start, start, [], new_seq))  

insertion_converter(sequence)

for i in insert_list:
    variant_dict.append(dict(zip(keys, i)))




deletions = []
occupied_deletion_range=set()
#deletion position generator
def mutater_deletion_pos(sequence_to_mutate_deletion):
    i=0
    while i < 10:
        deletion_length=random.randint(1,10)
        deletion_pos_start=random.randint(0, len(sequence_to_mutate_deletion)-1)
        deletion_pos_finish=min(deletion_pos_start + deletion_length, len(sequence_to_mutate_deletion)-1)
        
        deletion_range = range(deletion_pos_start, deletion_pos_finish)
        
        if any(position in snp_position for position in deletion_range):
            continue
        if any(position in deletion_range for position in occupied_insertion_range):
            continue
        if any(position in occupied_deletion_range for position in deletion_range):
            continue
        deletions.append((deletion_pos_start, deletion_pos_finish))
        occupied_deletion_range.update(deletion_range)
        i=i+1
    return deletions

mutater_deletion_pos(sequence)
#print(deletions)
        
        
                                       
#deletion mutater
deletion_list = []
def deletion_converter(sequence_to_add_deletion):
    for start, end in deletions:
        sequence_to_delete = []
        for i in range(start, end+1):
            base_delete = sequence_to_add_deletion[i]
            sequence_to_delete.append(base_delete)
        deletion_list.append(("deletion", start, end, sequence_to_delete, []))
deletion_converter(sequence)
#print(deletion_list)
        
for i in deletion_list:
    variant_dict.append(dict(zip(keys, i)))



#print(variant_dict[:5])



sorted_dictionary = sorted(variant_dict, key= lambda x: x["start"], reverse=True)
print(sorted_dictionary)

#apply the mutations in the dictionary to the genome
sequence_string = "".join(sequence)
mutated_sequence = sequence_string

for mutation in sorted_dictionary:
    mtype = mutation["mutation_type"]
    start = mutation["start"]
    end = mutation["end"]

    
    if mtype == "insert":
        new_seq = "".join(mutation["new"])
        mutated_sequence = (
            mutated_sequence[:start] +
            new_seq +
            mutated_sequence[start:]
        )

    
    elif mtype == "deletion":
        mutated_sequence = (
            mutated_sequence[:start] +
            mutated_sequence[end+1:]
        )

    
    elif mtype == "snp":
        new_base = "".join(mutation["new"])
        mutated_sequence = (
            mutated_sequence[:start] +
            new_base +
            mutated_sequence[start+1:]
        )
                            
    
#print(mutated_sequence[32223:32228+1])
#adding 1 to the index for samtools
variant_dict_for_samtools = copy.deepcopy(sorted_dictionary)
for entry in variant_dict_for_samtools:
    entry["start"] +=1
    entry["end"] +=1


with open("NC_037282_mutated.fasta", "w") as Fasta_file:
    Fasta_file.write(">Ecoli_mutated\n")
    for i in range(0, len(mutated_sequence), 60):
        Fasta_file.write(mutated_sequence[i:i+60] + "\n")
        
with open("NC_037282_mutant_dictionary.json", "w") as Dict_file:
    json.dump(variant_dict_for_samtools, Dict_file, indent=2)
    
    