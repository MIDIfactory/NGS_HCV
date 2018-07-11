#coding:utf-8
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import argparse
import pickle
from bwa_mapping import file_sam_to_be_processed
import numpy as np
import pandas as pd
import re
import itertools
import warnings
from datetime import datetime
from Bio import pairwise2
from blast_and_genotype_reads import get_genotype_of_the_reads_to_be_mapped
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
np.set_printoptions(threshold=np.nan)



############################## importo argomenti dal file principale #######################
scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) # + "/scripts"
dictionary_args = pickle.load(open( "%s/dictionary_args" % (scripts_folder),"r"))
################################################# funzioni ############################################################################
if sys.argv[0].split("/")[-1] == "processamento_SAM.py":
	print "______________________________________________________________"
	print "start of the sam processing..."

file_sam_to_be_processed=file_sam_to_be_processed()
# funzione per leggere il sam file come dataframe
def parse_sam():
	#selezioniamo il file sam che dobbiamo processare 
	filesam = open("./%s" % (file_sam_to_be_processed),"r")
	sam_raw=[x.strip().split("\t") for x in filesam] # each x is equal to a record(riga) of the sam file. lista sam non trattata,grezza.
	# trattiamo la lista delle records in modo che includano solo fino alla 11 colonna(in caso c'è ne siano di più il resto lo togliamo tanto è opzionale)
	sam=[]
	for record in sam_raw:
		sam.append(record[0:11])
	# i campi fino al phred score sono obbligatori ci devono essere in qualsiasi sam. Ci sono "11 mandatory field" nel sam
	column_names=["reads_id","sam_flag","ref_name","left_pos_map","mapping_quality","cigar","ref_name_R2","left_pos_map_R2","template_lenght","read_seq","phred_score"]
	samdf=pd.DataFrame(sam,columns=column_names) # trattiamo il sam come dataframe
	samdf=samdf.loc[2:] # togliamo le righe che riportano l'header,quindi consideriamo il dataframe a partire dalla terza riga(che ha label index pari a 2)
	samdf.reset_index(inplace=True,drop=True) 
	# di default reset index crea un nuovo dataframe se voglio che agisca sul vecchio dataframe uso inplace=True
	# drop=True evita che la vecchia serie di index mi sia aggiunta come colonna al dataframe
	return samdf

# significato colonne file sam:
# - sam_flag: indica se come si sono mappate le reads,se si sono mappate oppure no ecc... ogni numero ha un preciso significato
#	-> nel caso di pair-end mapping 
#	----83  indicano che quella read è la prima read del pair, si è mappata correttamente,tutta nell'insert size,nel reverse strand e la sua mate(la seconda read della coppia) si è mappata nel forward (strand) della ref . La sua mate è la read dopo 
#	---- 163 indica che quella read è la seconda read del pair,si è mappata correttamente,tutta nell'insert size,nel forward strand(strand +) mentre la sua mate si è mappata nel reverse strand(-). la sua mate è la read precedente  
#	---- 77 e 141 significa che entrambe le reads non hanno mappato
#	---- 81 e 161 stessa situazione di 83 e 163 ma sono uscite fuori dall'insert size
#	---- alcune volte potrebbero esserci dei flag strani tipo 2161 che indicano supplementary alignment,cioè che la read è stata spezzata in due,per esempio in caso di chromosomal inversion in e la read mappa su un punti di rottura allora un pezzo mappa in uno strand e uno nell'altro.
#	-> nel caso di single end mapping
#	----i flag 0 e 16 significa che la read si è mappata correttamente. 0 significa che la read si è mappata sullo strand forward della reference,quindi quello fornito. 16 significa che la read si è mappata sul reverse strand della reference.
#	---- 4 significa che la read non si è mappata
# - TLEN:TLEN(observed template lenght) is the distance from the 5'-most to 3'-most position.(quindi in pratica la lunghezza della read a cui si riferisce la riga. un valore negativo si riferisce al fatto che la read si è mappata al reverse strand,un valore positivo significa che è si è mappata sul forward strand della reference.


# funzione per ottenere informazioni varie sull'output del sam
def sam_reads_information(output): # con output stabilisco se portare fuori qualche informazione come variabile
	samdf=parse_sam()
	#---------------------- informazioni generali comuni sia con il single end che con il pair end mapping --------------------------
	# quante reads in totale sono state processate 
	n_tot_reads=samdf.shape[0]

	# quante reads non sono state mappate
	samdf_no_mapped=samdf[samdf["ref_name"]=="*"]
	n_reads_not_map=samdf_no_mapped.shape[0]

	# quante reads sono state mappate in generale
	samdf_mapped=samdf[(samdf["ref_name"]!="*")  & (samdf["cigar"]!="*") & (samdf["left_pos_map"]!="*") ]
	n_reads_map=samdf_mapped.shape[0]
	#---------------------- informazioni specifiche in base al tipo di mapping eseguito(flag del mapping),single end o pair end ----------------------
	# nel caso del pair-end mapping
	if dictionary_args["paired_end_mapping"]==True:
		# quante reads sono state entrambe mappate correttamente e dentro l'insert size
		flag_list=[99,147,83,163] # flag da tenere in considerazione
		samdf_mapped_correct=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_map_correct=samdf_mapped_correct.shape[0]
		# quante reads sono state mappate all'interno dell'insert size ma con orientamento sbagliato
		flag_list=[67,131,115,179] # flag da tenere in considerazione
		samdf_mapped_wrong_orientation=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_map_wrong_orientation=samdf_mapped_wrong_orientation.shape[0]
		# quante reads sono state mappate in maniera univoca,quindi con il giusto orientamento,ma fuori dall'insert size
		flag_list=[81,161,97,145,65,129,113,177] # flag da tenere in considerazione
		samdf_mapped_out_insert_size=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_map_out_insert_size=samdf_mapped_out_insert_size.shape[0]
		# quante reads sono mappate in cui solo un elemento della coppia è stata mappata
		flag_list=[73,133,89,121,165,181,101,117,153,185,69,137] # flag da tenere in considerazione
		samdf_mapped_1_map=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_map_1_map=samdf_mapped_1_map.shape[0]
		# quante reads sono state entrambe mappate correttamente, più quelle in cui solo una ha mappato correttamente,ma tutte cmq dentro l'insert size 
		flag_list=[99,147,83,163,73,133,89,121,165,181,101,117,153,185,69,137] # flag da tenere in considerazione
		samdf_mapped_correct_plus_1_map=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_mapped_correct_plus_1_map=samdf_mapped_correct_plus_1_map.shape[0]
		# quante reads sono state entrambe mappate correttamente,quelle in cui solo una ha mappato correttamente,e quelle mappate con orientamento sbagliato ma tutte cmq dentro l'insert size,  
		flag_list=[99,147,83,163,73,133,89,121,165,181,101,117,153,185,69,137,67,131,115,179] # flag da tenere in considerazione
		samdf_mapped_correct_all=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_mapped_correct_all=samdf_mapped_correct_all.shape[0]
		# quante reads sono state entrambe mappate correttamente,quelle in cui solo una ha mappato correttamente, quelle mappate con orientamento sbagliato,dentro l'insert size piu quelle che hanna mappato fuori dall'insert size  
		flag_list=[99,147,83,163,73,133,89,121,165,181,101,117,153,185,69,137,67,131,115,179,81,161,97,145,65,129,113,177] # flag da tenere in considerazione
		samdf_mapped_all_types=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_mapped_correct_all=samdf_mapped_all_types.shape[0]
		
		
	# nel caso del single end mapping
	elif dictionary_args["paired_end_mapping"]==False:
		# quante reads sono state mappate correttamente e dentro l'insert size. Nel caso del single end in realtà non c'è ne bisogno perchè i flag che ci interessa già vengono selezionati escludendo i "*" da ref_name	
		# Nel caso del single end in realtà non c'è ne bisogno perchè i flag che ci interessano già vengono selezionati escludendo i "*" da ref_name			
		flag_list=[0,16] # flag da tenere in considerazione
		samdf_mapped_all_types=samdf_mapped.loc[samdf_mapped['sam_flag'].astype(int).isin(flag_list)]
		n_reads_map_correct=samdf_mapped_all_types.shape[0] # numero di reads mappate correttamente
	# per il momento scegliamo come reads correttemente mappate tutti i tipi di reads mappate(se siamo nel pair_end,nel single end non c'è questa distinzione ovviamente)
	samdf_mapped_correct=samdf_mapped_all_types
	# -------------------- selezioniamo delle reads escludendo i cigar strani(presenza di N) -------------------------------------------
	cigar_list_correct=[]
	for cigar in samdf_mapped_correct["cigar"]:
		# escludiamo i cigar N,perchè,come scritto nella documentazione del SAM, rappresenta un introne 
		#e quindi può capitare solo se allinei m_rna con genomi eucariotici. Non ha senso in questo contesto e printa un warning tra l'altro
		if  "N" in cigar:
			warnings.warn("Found N cigar during the scan of the alignment. It indicates intron presence or unknown mapping,please check your input reads",RuntimeWarning)
			continue
		else:
			cigar_list_correct.append(cigar)
	# samdf_mapped_correct contiene tutte le informazioni sui mapping avvenuti correttamente(in caso del pair-end,entrambe mappate e dentro l'insert size) senza cigar strani.
	samdf_mapped_correct=samdf_mapped_correct.loc[samdf_mapped_correct['cigar'].astype(str).isin(cigar_list_correct)]
	n_reads_mapped_correct_and_analyzed=samdf_mapped_correct.shape[0]
	# prendiamo solo quelle mappate correttamente e il loro numero
	if output=="mapped_correct":
		return samdf_mapped_correct,n_reads_mapped_correct_and_analyzed
	# prendiamo tutte le informazioni
	elif output=="all":
		return n_tot_reads,n_reads_not_map,n_reads_map,n_reads_mapped_correct_and_analyzed
	# prendiamo il genotipo a cui appartengono le reads(ci serve dopo per selezionare l porzione del rule set corretta)
	elif output=="genotype_ref_map":
		if dictionary_args["user_ref"]==None:
			# se usano le ref del software il nome del genotipo possiamo anche prenderlo dal sam(nel nome della ref è indicato il genotipo scelto dal blast)
			# in questo caso non c'è bisogno di distinguere fra caso in lo user conosce il genotipo e caso in cui non lo conosce,perchè tanto l'informazione giusta sul genotipo sarà cmq contenuta nel nome della ref presente nel sam
			name_ref_of_mapping=samdf_mapped_correct["ref_name"].iloc[0]
			genotype_of_ref_of_mapping=name_ref_of_mapping.split("_")[1]
		elif dictionary_args["user_ref"]!=None:
			# se invece usiamo la ref dello user,possiamo prendere il genotipo solo dal blastout,perchp come nome ref nel sam c'è quella dello user che potrebbe non avere indicato il genotipo nel nome della ref.
			# in questo caso distinguiamo fra caso in cui il genotipo delle reads è sconosciuto allo user e caso in cui è conosciuto
			if dictionary_args["genotype"]=="unknown": # caso in cui è sconosciuto
				diz_genotype=get_genotype_of_the_reads_to_be_mapped(get_genotype=True,print_inf=False,remove_blastout=False,check_genotype=False)
				genotype_of_ref_of_mapping=diz_genotype["blast_genotype"] # è stato fatto il blast per sapere il genotipo
			elif dictionary_args["genotype"]!="unknown":
				genotype_of_ref_of_mapping=dictionary_args["genotype"]	
		return genotype_of_ref_of_mapping


# funzione per analizzare e controllare le informazioni ricavate dal mapping
def check_reads_information():
	n_tot_reads,n_reads_not_map,n_reads_map,n_reads_mapped_correct_and_analyzed=sam_reads_information(output="all")
	#------------ printiamo la percentuale di reads mappate correttamente e analizzate rispetto al totale delle reads
	reads_correct_mapped_percent=float(n_reads_mapped_correct_and_analyzed)/float(n_tot_reads) 
	print "---------------------------------------------------------"
	print "percentage of reads correctly mapped:%s" % (round(reads_correct_mapped_percent,3))
	if reads_correct_mapped_percent < 0.5 :
		warnings.warn("\n:::::::::::\nlow mapping rate!, less than 50 percent (%s mapped out of %s total reads)\n:::::::::::" % (n_reads_mapped_correct_and_analyzed,n_tot_reads),RuntimeWarning)	
	print "---------------------------------------------------------"
	# controlliamo che la percentuale di reads mappate non sia 0,in caso usciamo dall'esecuzione
	if n_reads_mapped_correct_and_analyzed==0:
		raise RuntimeError("No reads mapped,possible reasons:\n1)the input reads files has low quality(experiment issues)\n2)input gene provided by the user is wrong\n3) if a user_ref has been chosen, the user_ref may be inadapt")


# funzione per ottenere le reference tradotte e posizioni associate rispetto alla reference h77,ovvero il genoma del genotipo 1a.
# infatti le mutazioni riportate da geno2pheno si riferiscono agli aminoacidi.
def reference_selection_and_trad(translation):
	samdf_mapped_correct,_=sam_reads_information(output="mapped_correct")
	genotype_of_ref_of_mapping=sam_reads_information(output="genotype_ref_map")
	databaseHCV_folder=os.path.abspath("/".join((os.path.abspath("/".join(sys.argv[0].split("/")[:-1]))).split("/")[:-1])) + "/databases"
	# reference h77 is always present because i need it to insert the position of the mutation indicated in the rule set,that is based on h77 - c'è scritto nell'articolo
	# h77 si riferisce al genotipo 1a,ovviamente dovrò considerare la reference in base al gene NS5A,NS5B,NS3

	#---------------- caso in cui lo user sceglie di usare le ref dentro il software -------------
	if dictionary_args["user_ref"]==None:
		# -------------------- importiamo il multiple sequence alignment delle references tradotte o nucleotidiche -----------------------------------------------------
		# importiamo il multiple sequence alignment eseguito con clustal fra tutte le reference genotipiche a disposizione di un certo gene
		if translation==True:
			SeqIO_input="%s/databaseHCV_bwa/Alignment_references_all_genotypes/HCV_alignment_reference_%s_all_genotypes_protein.fst" % (databaseHCV_folder,dictionary_args["gene"]) 
		elif translation==False:
			SeqIO_input="%s/databaseHCV_bwa/Alignment_references_all_genotypes/HCV_alignment_reference_%s_all_genotypes.fst" % (databaseHCV_folder,dictionary_args["gene"]) 
		alignment_references_all_genotypes=SeqIO.parse(SeqIO_input,"fasta")
		alignment_references_all_genotypes=pd.DataFrame(alignment_references_all_genotypes,index=["1a","1b","2","3","4"])
		#come già detto la reference del genotipo 1a,ovvero h77, è quella a cui si riferiscono le posizioni del rule set di geno2pheno.
		#seleziono inoltre, insieme alla ref h77,la reference contro cui le reads sono state mappate
		ref_1a_plus_ref_of_mapping=alignment_references_all_genotypes.loc[["1a","%s" % (genotype_of_ref_of_mapping)]].astype(str).T
		if ref_1a_plus_ref_of_mapping.columns.values[0] != ref_1a_plus_ref_of_mapping.columns.values[1]:	 
			# in pratica con il codice di sotto elimino dal dataframe i casi in cui ho "-" sia nella ref h77 sia nel ref di mapping perchè non ha senso che ci siano,ci sono solo perchè ho fatto un multiple alignment quindi quei "-" si riferiscono a terze sequenze.
			ref_1a_plus_ref_of_mapping=ref_1a_plus_ref_of_mapping.loc[(ref_1a_plus_ref_of_mapping['1a'] != "-") | (ref_1a_plus_ref_of_mapping['%s' % (genotype_of_ref_of_mapping)] != "-")]
			ref_1a_plus_ref_of_mapping.reset_index(inplace=True,drop=True)
			ref_1a=ref_1a_plus_ref_of_mapping["1a"]
			ref_of_mapping=ref_1a_plus_ref_of_mapping["%s" % (genotype_of_ref_of_mapping)]
			array_index_ref_1a=ref_1a.index.values
		# nel caso le reads abbiano genotipo 1a,avremmo un dataframe con colonne con lo stesso nome e questo non va bene.
		elif ref_1a_plus_ref_of_mapping.columns.values[0] == ref_1a_plus_ref_of_mapping.columns.values[1]:
			new_columns_name=ref_1a_plus_ref_of_mapping.columns.values
			new_columns_name[1]="1a_ref_of_mapping"
			ref_1a_plus_ref_of_mapping.columns=new_columns_name
			# in pratica con il codice di sotto elimino dal dataframe i casi in cui ho "-" sia nella ref h77 sia nel ref di mapping perchè non ha senso che ci siano,ci sono solo perchè ho fatto un multiple alignment quindi quei "-" si riferiscono a terze sequenze.
			ref_1a_plus_ref_of_mapping=ref_1a_plus_ref_of_mapping.loc[(ref_1a_plus_ref_of_mapping['1a'] != "-") | (ref_1a_plus_ref_of_mapping['1a_ref_of_mapping'] != "-")]
			ref_1a_plus_ref_of_mapping.reset_index(inplace=True,drop=True)
			ref_1a=ref_1a_plus_ref_of_mapping["1a"]
			ref_of_mapping=ref_1a_plus_ref_of_mapping["1a_ref_of_mapping"]
			array_index_ref_1a=ref_1a.index.values
	# ----------------------- caso in cui lo user decide di usare la sua reference -------------------------------------------
	elif dictionary_args["user_ref"]!=None:
		SeqIO_input_ref_1a="%s/databaseHCV_bwa/databaseHCV_reference_1a/HCV_1a_%s.TXT" % (databaseHCV_folder,dictionary_args["gene"]) 
		ref_1a=SeqIO.read(SeqIO_input_ref_1a,"fasta")
		user_ref=SeqIO.read("%s" % (dictionary_args["user_ref"]),"fasta")
		if translation==False:
			needle_cline = NeedleCommandline(asequence=SeqIO_input_ref_1a, bsequence="%s" % (dictionary_args["user_ref"]),gapopen=10, gapextend=0.5,outfile="needle.txt")
			stdout,stderr = needle_cline()
		elif translation==True:
			with warnings.catch_warnings():
				warnings.simplefilter("ignore") # ignoriamo l'inutile biopython warning sul trimming codons.
				ref_1a=str(ref_1a.seq.translate())			 
				user_ref=str(user_ref.seq.translate())
			with open("ref_1a_needle_input.fasta", "w") as ref_1a_needle_input_file:
				ref_1a_needle_input_file.write(">ref_1a\n%s" % (ref_1a))
			with open("user_ref_needle_input.fasta", "w") as user_ref_needle_input_file:
				user_ref_needle_input_file.write(">user_ref\n%s" % (user_ref))
			needle_cline = NeedleCommandline(asequence="./ref_1a_needle_input.fasta", bsequence="./user_ref_needle_input.fasta",gapopen=10, gapextend=0.5,outfile="needle.txt")
			stdout,stderr = needle_cline()
			# rimuoviamo i needle input files appena creati
			os_input="rm -f ./ref_1a_needle_input.fasta;rm -f ./user_ref_needle_input.fasta"
			os.system(os_input)

		alignment_ref_1a_vs_user_ref=AlignIO.read("needle.txt", "emboss")
		ref_1a_seq_aln=np.array(list(alignment_ref_1a_vs_user_ref[0].seq))
		user_ref_seq_aln=np.array(list(alignment_ref_1a_vs_user_ref[1].seq))
		
		df_ref_1a_seq_aln=pd.Series(ref_1a_seq_aln)
		df_ref_1a_seq_aln=df_ref_1a_seq_aln.rename("1a_ref")
		df_user_ref_seq_aln=pd.Series(user_ref_seq_aln)
		df_user_ref_seq_aln=df_user_ref_seq_aln.rename("user_ref")
		ref_1a_plus_ref_of_mapping=pd.concat((df_ref_1a_seq_aln,df_user_ref_seq_aln),axis=1)
		# in pratica con il codice di sotto elimino dal dataframe i casi in cui ho "-" sia nella ref h77 sia nel ref di mapping perchè non ha senso che ci siano,ci sono solo perchè ho fatto un multiple alignment quindi quei "-" si riferiscono a terze sequenze.
		ref_1a_plus_ref_of_mapping=ref_1a_plus_ref_of_mapping.loc[(ref_1a_plus_ref_of_mapping['1a_ref'] != "-") | (ref_1a_plus_ref_of_mapping['user_ref'] != "-")]

		array_index_ref_1a=df_ref_1a_seq_aln.index.values
		ref_1a=ref_1a_plus_ref_of_mapping["1a_ref"]
		ref_of_mapping=ref_1a_plus_ref_of_mapping["user_ref"]
		array_index_ref_1a=ref_1a.index.values
		#---- checkpoint !!! -----------
		# assicuriamoci che il trattamento della ref dello user sia andato a buon fine,la dimensione dei due array deve essere uguale
		if ref_1a.shape != ref_of_mapping.shape:
			raise RuntimeError("ref_1a and user_ref_of_mapping have different shape,check the code")
		#-------------------------------------
	
	# ----------------------- posizione e aminoacido o nucleotide associato della read-----------------------------
	# le posizioni di reference diverse dalla 1a sono calibrate in maniera che corrispondono alle posizioni di 1a tradotta 
	# perchè le posizioni delle mutazioni aminoacidiche di geno2pheno si riferisono alle posizioni di 1a,ovvero h77.
	new_list_index_ref_of_mapping=[]
	new_list_index_ref_1a=[]
	insertion_index_ref_mapping=0
	insertion_index_ref_1a=0
	for index in array_index_ref_1a:
		ref_1a_aa=ref_1a[index]
		ref_of_mapping_aa=ref_of_mapping[index]
		# caso in cui non c'è inserzione nè delezione fra la ref h77,cioè 1a  e la ref su cui abbiamo mappato le reads
		#l'index della reference di mappaggio è uguale all'index di h77
		if ((ref_1a_aa!="-" and ref_of_mapping_aa!="-") or (ref_1a_aa=="-" and ref_of_mapping_aa=="-")):
			if insertion_index_ref_mapping==0: 
				new_list_index_ref_of_mapping.append("%s" % str(index+1)) # inserisco +1 perchè voglio che le posizioni partano da 1 seguendo il modo di indexare di del mapping
				new_list_index_ref_1a.append("%s" % str(index+1))
			else:
				new_list_index_ref_of_mapping.append("%s" % str(index+1-insertion_index_ref_mapping))
				new_list_index_ref_1a.append("%s" % str(index+1-insertion_index_ref_1a))

		# caso in cui c'è una delezione nella ref di mappaggio rispetto alla ref h77. 
		# Significa che nella ref h77 c'è un inserzione rispetto alla ref di mappaggio
		# l'index della reference di mappaggio è diverso dall'index di h77 
		elif ref_1a_aa!="-" and ref_of_mapping_aa=="-":
			insertion_index_ref_mapping+=1
			new_list_index_ref_of_mapping.append("%s+%s" % (str(index+1-insertion_index_ref_mapping),str(insertion_index_ref_mapping)))
			if insertion_index_ref_1a==0:
				new_list_index_ref_1a.append("%s" % str(index+1))
			elif insertion_index_ref_1a!=0:
				new_list_index_ref_1a.append("%s" % str(index+1-insertion_index_ref_1a))			
		# caso in cui c'è un inserzione nella ref di mappaggio rispetto a h77. 
		# Significa che nella ref h77 c'è una delezione.
		# l'index della reference di mappaggio è uguale rispetto all'index della ref h77
 		elif ref_1a_aa=="-" and ref_of_mapping_aa!="-":
			insertion_index_ref_1a+=1
			new_list_index_ref_1a.append("%s+%s" % (str(index+1-insertion_index_ref_1a),str(insertion_index_ref_1a)))
			if insertion_index_ref_mapping==0: 
				new_list_index_ref_of_mapping.append("%s" % str(index+1))
			elif insertion_index_ref_mapping!=0:
				new_list_index_ref_of_mapping.append("%s" % str(index+1-insertion_index_ref_mapping))
					
	array_index_ref_of_mapping=np.array([new_list_index_ref_of_mapping])
	array_aa_ref_of_mapping=np.array([ref_of_mapping])
	array_index_vs_aa_ref_of_mapping=np.concatenate((array_index_ref_of_mapping,array_aa_ref_of_mapping),axis=0).T
	array_index_ref_1a=np.array([new_list_index_ref_1a]) # se voglio che l'array_index_ref_1a contenga gli index senza considerare le inserzioni della ref elimino questo codice.
	# gli index però sono ancora quelli della ref_of_mapping. ci portiamo dietro le posizioni della ref_1a nell'array_index_ref_1a
	return array_index_vs_aa_ref_of_mapping,array_index_ref_1a 



# in pratica array_index_vs_aa_ref_of_mapping contiene come prima colonna le posizioni aminoadiche riferite alla ref h77(ref 1a) della ref di mapping
# quando andremo ad analizzare le reads mappate,compareremo le posizioni delle reads e aminoacidi associati con le posizione e aa riportati in questo array

# funzione per definire il frame quindi trovare la posizione in frame iniziale della read
def frame_definition():
	samdf_mapped_correct,_=sam_reads_information(output="mapped_correct") # prendiamo le reads correttamente mappate
	array_left_pos_map=np.array(samdf_mapped_correct["left_pos_map"]).astype(int)
	array_pos_start_read_in_frame_nucleot=np.array([]) # array dove riporto le correzioni da effetuare in base al frame  prima di tradurre
	array_pos_start_read_in_frame_aa=np.array([]) # array dove riporto le correzioni da effettuare quando devo ricarvarmi le posizioni delle seq già tradotte
	for pos in array_left_pos_map:
		# caso in cui la read è gia in frame. Se l'inizio mappato della read coincide con l'inizio della tripletta,allora la left pos map della read sarà un numero divisibile per 3,
		if (pos+2) % 3 == 0 :
			# frame=1
			pos_start_read_in_frame=0
			array_pos_start_read_in_frame_nucleot=np.append(array_pos_start_read_in_frame_nucleot,pos_start_read_in_frame)
			array_pos_start_read_in_frame_aa=np.append(array_pos_start_read_in_frame_aa,2) # perchè se vuoi avere la corrispondente posizione aminoacidica delle posizioni nucleotidiche devi andare alla fine della tripletta.
		# NOTA: non ho bisogno in realtà in questo caso di queste variabili perchè io già parto a tradurre dalla sequenza nucleotidica già in frame,quindi devo sempre spostarmi di due.
		# gli altri due casi in cui la read non è subito in frame
		elif (pos+1) % 3 == 0 :
			# frame=2
 			pos_start_read_in_frame=2
			array_pos_start_read_in_frame_nucleot=np.append(array_pos_start_read_in_frame_nucleot,pos_start_read_in_frame)
			array_pos_start_read_in_frame_aa=np.append(array_pos_start_read_in_frame_aa,4)

		elif (pos) % 3 == 0 :
			# frame=3
 			pos_start_read_in_frame=1
			array_pos_start_read_in_frame_nucleot=np.append(array_pos_start_read_in_frame_nucleot,pos_start_read_in_frame)
			array_pos_start_read_in_frame_aa=np.append(array_pos_start_read_in_frame_aa,3)
	return array_pos_start_read_in_frame_nucleot,array_left_pos_map,array_pos_start_read_in_frame_aa



def extract_cigar_informations():
	samdf_mapped_correct,_=sam_reads_information(output="mapped_correct") # prendiamo le reads correttamente mappate	
	array_cigar=np.array(samdf_mapped_correct["cigar"]).astype(str)
	#array_cigar=np.array(["5S46M1D3I25S","3S4D52M3I12D5S","3S30D66I35M7S","6S3D35M4D3I32M4D3I32M4D15S"]) # array inventato che uso come test ho fatto che le lunghezza del cigar corrisponda alla lunghezza totale delle prime 4 read del campione 30, combinazione batch 100 1 0
	lista_df_letters_vs_digits_all_reads=[] # lista dove riporto tutte i dataframes dei cigar. Lettere vs digits
	lista_string_cigar_informations_no_S_all_reads=[] # lista dove riporto tutte le lettere del cigar ripetute per il loro valore di digit
	lista_string_cigar_informations=[]
	for cigar in array_cigar:
		array_letters=np.array([list(itertools.ifilterfalse(str.isdigit,cigar))])
		array_digits=np.array([re.findall('\d+', cigar)])
		array_letters_vs_digits=np.concatenate((array_letters,array_digits),axis=0)
		df_letters_vs_digits=pd.DataFrame(array_letters_vs_digits,columns=array_letters_vs_digits[0])
		df_letters_vs_digits.drop(df_letters_vs_digits.index[0],inplace=True)
		df_letters_vs_digits=df_letters_vs_digits.astype(int)
		lista_df_letters_vs_digits_all_reads.append(df_letters_vs_digits)
		string_cigar_informations=""
		string_cigar_informations_no_S=""
		for i in range(len(df_letters_vs_digits.columns.values)):
			value=df_letters_vs_digits.iloc[:,i][1]
			if df_letters_vs_digits.columns.values[i]=="S":
				string_S="S"*value
				string_cigar_informations=string_cigar_informations+string_S
			elif df_letters_vs_digits.columns.values[i]=="M":
				string_M="M"*value
				string_cigar_informations=string_cigar_informations+string_M
				string_cigar_informations_no_S=string_cigar_informations_no_S+string_M
			elif df_letters_vs_digits.columns.values[i]=="D":
				string_D="D"*value
				string_cigar_informations=string_cigar_informations+string_D
				string_cigar_informations_no_S=string_cigar_informations_no_S+string_D
			elif df_letters_vs_digits.columns.values[i]=="I":
				string_I="I"*value
				string_cigar_informations=string_cigar_informations+string_I
				string_cigar_informations_no_S=string_cigar_informations_no_S+string_I
			elif df_letters_vs_digits.columns.values[i]=="=":
				string_eg="="*value
				string_cigar_informations=string_cigar_informations+string_eg
				string_cigar_informations_no_S=string_cigar_informations_no_S+string_eg
			elif df_letters_vs_digits.columns.values[i]=="X":
				string_X="X"*value
				string_cigar_informations=string_cigar_informations+string_X
				string_cigar_informations_no_S=string_cigar_informations_no_S+string_X

		# se incontra H non cambia nulla se considerarla o meno perchè tanto non consuma di posizioni di read(dove viene già tolta in automatico) nè posizioni di ref. Vedi esempio sam pdf
		# possiamo anche ignorare P,perchè P è una mutazione silente quindi non ha senso considerarla(in ogni caso non consuma ne reference ne read)
		lista_string_cigar_informations_no_S_all_reads.append(string_cigar_informations_no_S)
		lista_string_cigar_informations.append(string_cigar_informations)
	return lista_df_letters_vs_digits_all_reads,lista_string_cigar_informations_no_S_all_reads,lista_string_cigar_informations
			


#concettualmente:
# reference_portion: dove c'è M e D  ci sono lettere nella reference,in corrispondenza di D c'è trattino nella read
# read_portion: dove c'è M e I ci sono lettere nella read, in corrispondenza di I ci sono trattini nella reference

def extract_ref_portion_mapped():
	dataframe_ref_portion_list=[] # lista di dataframe dove metto tutti i dataframe delle ref portion mappate dalle reads
	lista_df_letters_vs_digits_all_reads,lista_string_cigar_informations_no_S_all_reads,_=extract_cigar_informations()
	array_pos_start_read_in_frame_nucleot,array_left_pos_map,array_pos_start_read_in_frame_aa=frame_definition()
	array_pos_start_read_in_frame_nucleot=array_pos_start_read_in_frame_nucleot.astype(int)
	array_index_vs_aa_ref_of_mapping,_=reference_selection_and_trad(False)
	df_index_vs_aa_ref_of_mapping=pd.DataFrame(array_index_vs_aa_ref_of_mapping,columns=["pos","nucl"])
	df_index_vs_aa_ref_of_mapping.set_index("pos",inplace=True)
	for i in range(len(lista_string_cigar_informations_no_S_all_reads)):
		string_cig_x=lista_string_cigar_informations_no_S_all_reads[i]
		pos_start_ref_portion_mapped=array_left_pos_map[i]+array_pos_start_read_in_frame_nucleot[i] # array_left_pos_map[-1]  è sbagliato ovviamente ma questo è una prova
		array_pos=np.array([])
		l=-1 # conta le lettere
		ins=0 # conta le inserzioni
		for letter in string_cig_x:
			l+=1
			if letter!="I":
				if ins==0:
					array_pos=np.append(array_pos,str(int(l)))
				elif ins!=0:
					array_pos=np.append(array_pos,str(int(l-ins)))	
			elif letter=="I":
				ins+=1
				array_pos=np.append(array_pos,"%s+%s" % (l-ins,ins))

		#------------------------------------------ sistemiamo il frame delle pos della ref --------------------------------------

		# In questo caso ho già preso le string cigar information senza la S, non ho bisogno di sottolineare che deve saltare le S tanto non ci sono.
		array_pos_of_ref_mapped=np.array([])
		for pos in array_pos:
			if "+" not in pos:
				pos_int=int(float(pos)) # non serve aggiungere il +1 perchè anche se parte a contare da 0,noi aggiungiamo la posizione di inizio del mapping,non usiamo la funzione range che non include la fine. Facendo l'addizione includiamo la fine.
				# cioè aggiungendo per esempio 1258,partiamo da 1258 non partiamo da 1257,perchè fare 0 + 1258 risulta 1258,non escludiamo l'ultimo numero,come sarebbe successo se avessimo usato range
				pos_int=pos_int+array_left_pos_map[i]
				array_pos_of_ref_mapped=np.append(array_pos_of_ref_mapped,str(int(pos_int)))
			elif "+" in pos:
				pos_int_I=int(pos.split("+")[0])
				pos_int_I=str(int(pos_int_I+array_left_pos_map[i]))
				addition=pos.split("+")[1]
				pos_int_I="+".join((pos_int_I,addition))
				array_pos_of_ref_mapped=np.append(array_pos_of_ref_mapped,pos_int_I)
		# sistemiamo il frame iniziale
		array_pos_of_ref_mapped=array_pos_of_ref_mapped[array_pos_start_read_in_frame_nucleot[i]:]
		# sistemiamo il frame finale
		lenght_tot_ref_mapped=array_pos_of_ref_mapped.shape[0]
		# lenght_tot_ref_mapped deve essere divisibile per tre,in caso contrario dobbiamo sistemare anche la posizione finale del porzione affinchè sia tutto in frame
		if (lenght_tot_ref_mapped-1) % 3 == 0:
			array_pos_of_ref_mapped=array_pos_of_ref_mapped[:lenght_tot_ref_mapped-1] # con la lunghezza calcolata con .shape arrivo fino alla fine per arrivare al nucleotide prima della fine sottraggo meno 1.
		elif (lenght_tot_ref_mapped-2) % 3 == 0:
			array_pos_of_ref_mapped=array_pos_of_ref_mapped[:lenght_tot_ref_mapped-2]	
		elif (lenght_tot_ref_mapped) % 3 == 0:
			array_pos_of_ref_mapped=array_pos_of_ref_mapped

		#--------------- mettiamo in dataframe --------------------------------
	
		array_pos_of_ref_mapped=np.array([array_pos_of_ref_mapped]).T
		array_ref_nucl_mapped=np.zeros((array_pos_of_ref_mapped.shape[0],1)).astype(str)
		array_pos_vs_nucl_ref_mapped=np.concatenate((array_pos_of_ref_mapped,array_ref_nucl_mapped),axis=1)
		df_pos_of_ref_mapped_vs_empy_nucl=pd.DataFrame(array_pos_vs_nucl_ref_mapped,columns=["pos","nucl"])
		df_pos_of_ref_mapped_vs_empy_nucl.set_index("pos",inplace=True)
		pos_values=df_pos_of_ref_mapped_vs_empy_nucl.index.values
		df_pos_of_ref_mapped_vs_nucl=df_index_vs_aa_ref_of_mapping.reindex(pos_values,fill_value="-") # si può usare anche il .loc ma il futureWarning mi ha avvertito che in futuro in presenza di Nan darà key error quindi mi ha suggerito di sostuirlo con reindex()
		dataframe_ref_portion_list.append(df_pos_of_ref_mapped_vs_nucl)
 
	return dataframe_ref_portion_list
	

# funzione per selezionare la parte di read da tenere in considerazione quando paragoniamo gli aminoacidi a quelli della porzione di ref su cui ha mappato
# bisogna togliere le basi soft clipping(S) dalla sequenza della read perchè nella sequenza dell'output del mapping appaiono.
# invece la lef_pos_map calcolata da bwa si riferisce sempre all'inizio di basi allineate già senza considerare il soft clipping.
def extract_read_portion_mapped():
	dataframes_ref_portion_vs_read_portion_list=[] # metteremo alla fine tutti i dataframes che comprano le read portion con la ref portion per ogni allineamento
	samdf_mapped_correct,_=sam_reads_information(output="mapped_correct") # prendiamo le reads correttamente mappate
	lista_df_letters_vs_digits_all_reads,lista_string_cigar_informations_no_S_all_reads,lista_string_cigar_informations=extract_cigar_informations()
	array_pos_start_read_in_frame_nucleot,array_left_pos_map,array_pos_start_read_in_frame_aa=frame_definition()
	array_pos_start_read_in_frame_nucleot=array_pos_start_read_in_frame_nucleot.astype(int)
	dataframe_ref_portion_list=extract_ref_portion_mapped()
	all_reads_seq=samdf_mapped_correct["read_seq"]
	# è la stessa cosa identica fatta per la ref_portion mapped tranne che qua mettiamo i trattini in corrispondenza delle D non delle I
	for i in range(len(lista_string_cigar_informations_no_S_all_reads)):
		string_cig_x=lista_string_cigar_informations_no_S_all_reads[i]
		pos_start_read_portion_mapped=array_left_pos_map[i]+array_pos_start_read_in_frame_nucleot[i] # ovviamente è uguale ala pos_start_read_portion_mapped
		array_pos=np.array([])
		l=-1 # conta le lettere
		d=0 # conta le delezioni
		for letter in string_cig_x:
			l+=1
			if letter!="D":
				if d==0:
					array_pos=np.append(array_pos,str(int(l)))
				elif d!=0:
					array_pos=np.append(array_pos,str(int(l-d)))	
			elif letter=="D":
				d+=1
				array_pos=np.append(array_pos,"%s+%s" % (l-d,d))

		#-------------------- sistemiamo il frame iniziale e finale delle pos_read -------------------------
		array_pos_of_read_mapped=np.array([])
		for pos in array_pos:
			if "+" not in pos:
				pos_int=int(float(pos)) # non serve aggiungere il +1 perchè anche se parte a contare da 0,noi aggiungiamo la posizione di inizio del mapping,non usiamo la funzione range che non include la fine. Facendo l'addizione includiamo la fine.
				# cioè aggiungendo per esempio 1258,partiamo da 1258 non partiamo da 1257,perchè fare 0 + 1258 risulta 1258,non escludiamo l'ultimo numero,come sarebbe successo se avessimo usato range
				pos_int=pos_int+array_left_pos_map[i]
				array_pos_of_read_mapped=np.append(array_pos_of_read_mapped,str(int(pos_int)))
			elif "+" in pos:
				pos_int_D=int(pos.split("+")[0])
				pos_int_D=str(int(pos_int_D+array_left_pos_map[i]))
				addition=pos.split("+")[1]
				pos_int_D="+".join((pos_int_D,addition))
				array_pos_of_read_mapped=np.append(array_pos_of_read_mapped,pos_int_D)
		# sistemiamo dunque le posizioni per il frame iniziale
		array_pos_of_read_mapped=array_pos_of_read_mapped[array_pos_start_read_in_frame_nucleot[i]:]
		lenght_tot_read_mapped=array_pos_of_read_mapped.shape[0]
		# lenght_tot_read_mapped deve essere divisibile per tre,in caso contrario dobbiamo sistemare anche la posizione finale del porzione affinchè sia tutto in frame
		if (lenght_tot_read_mapped-1) % 3 == 0:
			array_pos_of_read_mapped=array_pos_of_read_mapped[:lenght_tot_read_mapped-1] # con la lunghezza calcolata con .shape arrivo fino alla fine per arrivare al nucleotide prima della fine sottraggo meno 1.
		elif (lenght_tot_read_mapped-2) % 3 == 0:
			array_pos_of_read_mapped=array_pos_of_read_mapped[:lenght_tot_read_mapped-2]	
		elif (lenght_tot_read_mapped) % 3 == 0:
			array_pos_of_read_mapped=array_pos_of_read_mapped

			 
		array_pos_of_read_mapped=np.array([array_pos_of_read_mapped]).T
		array_read_nucl_mapped=np.zeros((array_pos_of_read_mapped.shape[0],1)).astype(str)
		array_pos_vs_nucl_read_mapped=np.concatenate((array_pos_of_read_mapped,array_read_nucl_mapped),axis=1)
		df_pos_of_read_mapped_vs_empy_nucl=pd.DataFrame(array_pos_vs_nucl_read_mapped,columns=["pos","nucl"])
		df_pos_of_read_mapped_vs_empy_nucl.set_index("pos",inplace=True)
		pos_values=df_pos_of_read_mapped_vs_empy_nucl.index.values # non mi interessa che le posizioni  siano alterate,le posizioni mi servivano solo per mettere i trattini al posto giusto
		# ------------- adesso processiamo la sequenza della read ----------------
		array_read_seq_x=np.array(list(all_reads_seq.iloc[i]))
		array_string_cigar_complete_x=np.array(list(lista_string_cigar_informations[i]))
		array_reads_seq_x_complete=np.array([]) # creaimo la read_seq che include i trattini

		l=-1
		d=0
		for letter in array_string_cigar_complete_x: # array_string_cigar_complete_x, iteriamo su questo che è la sequenza più lunga perchè include le D al suo interno
			l+=1
			if letter!="D":
				if d==0:
					array_reads_seq_x_complete=np.append(array_reads_seq_x_complete,array_read_seq_x[l])
				elif d!=0:
					array_reads_seq_x_complete=np.append(array_reads_seq_x_complete,array_read_seq_x[l-d])
			elif letter=="D":
				d+=1
				array_reads_seq_x_complete=np.append(array_reads_seq_x_complete,"-")
		#------------ checkpoint !!!--------------
		# la somma di M,S,I,=,X deve essere uguale alla lunghezza della sequenza della read senza trattini(non sono riportati nella sequenza del sam)
		# quindi la somma di M,S,I,=,X,D deve essere uguale alla lunghezza della sequenza che include i trattini
		if array_string_cigar_complete_x.shape != array_reads_seq_x_complete.shape:
			print "problems_in_reads_of_iteration:%d,the shapes %s %s don't match" % (i,array_string_cigar_complete_x.shape,array_reads_seq_x_complete.shape)
			raise Exception("RuntimeError:sum of M,S,I,=,X,D not match lenght SEQ with -")
		#------------------------------------------ 
		array_reads_seq_x_complete=np.array([array_reads_seq_x_complete])
		array_string_cigar_complete_x=np.array([array_string_cigar_complete_x])
		# togliamogli l's clipping
		array_cigar_plus_read_seq_x=np.concatenate((array_reads_seq_x_complete,array_string_cigar_complete_x),axis=0).T
		array_read_seq_x_no_S=array_cigar_plus_read_seq_x[array_cigar_plus_read_seq_x[:,1]!="S",0]
		#-------------------------------------- mettiamo in frame i nucleotidi della read -------------------------------				
		# combinazione errore: 1000 10 8
		# ricordiamoci anche di mettere in frame i nucleotidi,abbiamo già messo in frame le posizioni,ma dobbiamo farlo anche per i nucleotidi
		# nel caso della ref avevamo già l'assoziaone posizioni con nucl,quindi è basta usare reindex qua invece sistemiamo la cosa esplicitamente
		array_read_seq_x_no_S=array_read_seq_x_no_S[array_pos_start_read_in_frame_nucleot[i]:]

		# sistemiamo i nucleotidi tenendo in considerazioni le posizioni delle D, che conosciamo,perchè le abbiamo trovato prima,sono conservate in df_pos_of_read_mapped_vs_empy_nucl,le D sono nelle posizioni con il +.
		# le pos values sono messe in frame anche alla fine,infatti togliamo i nucleotidi di array_read_seq_x_no_S che non sono in frame alla fine
		array_read_seq_x_final=np.array([])
		p=-1	
		
		for pos in pos_values:
			p+=1
			if "+" not in pos:
				array_read_seq_x_final=np.append(array_read_seq_x_final,array_read_seq_x_no_S[p])
			elif "+" in pos:
				array_read_seq_x_final=np.append(array_read_seq_x_final,"-")
		
 	
		# adesso mettiamo la porzione di read con la porzione di ref nello state dataframe,avente come index le posizioni della ref.
		array_read_seq_x_final=np.array([array_read_seq_x_final])
		pos_ref_portion=np.array([dataframe_ref_portion_list[i].index.values]).T
		nucl_ref_portion=np.array([dataframe_ref_portion_list[i]["nucl"]])

		array_nucl_ref_nucl_read=np.concatenate((nucl_ref_portion,array_read_seq_x_final),axis=0).T
		array_pos_ref_nucl_ref_nucl_read=np.concatenate((pos_ref_portion,array_nucl_ref_nucl_read),axis=1)
		df_pos_ref_nucl_ref_nucl_read=pd.DataFrame(array_pos_ref_nucl_ref_nucl_read,columns=["pos_ref","nucl_ref","nucl_read"])
		df_pos_ref_nucl_ref_nucl_read.set_index("pos_ref",inplace=True)
		
	
		#---------- checkpoint !!!!-----

		# assicuriamoci che la lunghezza delle read e la lunghezza della ref sia divisibile per 3. così possiamo proseguire con la traduzione
		assert (df_pos_ref_nucl_ref_nucl_read.shape[0]) % 3 == 0, "RuntimeWarning: iteration %s,nucleotide lenght of the read and ref is not divisible for three,translation may be effect by errors." % (i)

		# --------------------------------------
		dataframes_ref_portion_vs_read_portion_list.append(df_pos_ref_nucl_ref_nucl_read)
		

	return dataframes_ref_portion_vs_read_portion_list


	

	

# funzione per tradurre la read nella porzione di ref su cui ha mappato 
def translation_of_the_ref_portion_mapped():
	lista_dataframes_pos_ref_read_trad=[] # lista di dataframes in cui riporto tutti i dataframes che si riferiscono ogni allineamento
	dataframes_ref_portion_vs_read_portion_list=extract_read_portion_mapped()
	for x in range(len(dataframes_ref_portion_vs_read_portion_list)):
		read_seq_x=dataframes_ref_portion_vs_read_portion_list[x]["nucl_read"]
		ref_portion_seq_x=dataframes_ref_portion_vs_read_portion_list[x]["nucl_ref"]
		# -------------- calcoliamo le posizioni aminoacidiche corrispodenti delle pos nucleotidiche -----------------------
		tripletta_ref_nucl=[]
		pos_ref_aa_seq_x=[]
		n=0
		p=-1
		ins=0
		for nucl in ref_portion_seq_x:
			n+=1
			tripletta_ref_nucl.append(nucl)
			if len(tripletta_ref_nucl)!=3:
				continue
			elif len(tripletta_ref_nucl)==3:
				p+=1 # ogni 3 nucleotidi abbiamo una posizione aminoacidica
				if "-" not in tripletta_ref_nucl:
					if ins==0:
						pos_ref_aa_seq_x.append(str(p))
					elif ins!=0:
						pos_ref_aa_seq_x.append(str(p-ins))
				elif "-" in tripletta_ref_nucl:
					ins+=1
					pos_ref_aa_seq_x.append("%s+%s" % (str(p-ins),ins))
				del tripletta_ref_nucl[:]
		
		# calcoliamoci la left_pos_map aminoacidica
		m=0
		for pos in ref_portion_seq_x.index.values:
			if "+" not in pos: # la prima posizione divisibile per tre che incontra è quella che mi dà la posizione aminoacidica di inizio.
				if (int(pos))%3==0:
					m+=1
					left_pos_map_aa=int(pos)/3
					if m==1:
						break

		correct_pos_ref_aa_seq_x=[]
		# adesso convertiamole in posizioni reali del mapping
		for pos in pos_ref_aa_seq_x:
			if "+" not in pos:
				correct_pos_ref_aa_seq_x.append(int(pos)+left_pos_map_aa)
			elif "+" in pos:
				correct_pos_ref_aa_seq_x.append("%s+%s" % (str(int(pos.split("+")[0])+left_pos_map_aa),pos.split("+")[1]))
			
		#---------------- troviamo gli aminoacidi della read --------------------------------------------
		nucleotides_read_seq_x=read_seq_x.values
		n=0 # conta i nucleotidi
		read_seq_x_aa=np.array([])
		tripletta_nucl=[] # lista dove metti i nucleotidi che costituiscono una tripletta
		for nucl in nucleotides_read_seq_x:
			n+=1
			tripletta_nucl.append(nucl)
			if len(tripletta_nucl)!=3:
				continue
			elif len(tripletta_nucl)==3:
				tripletta_joined="".join(tripletta_nucl)
				#if x==9:
					#print tripletta_joined,n
				tripletta_joined=Seq(tripletta_joined)
				if "-" not in tripletta_joined:
					aa=tripletta_joined.translate()
					read_seq_x_aa=np.append(read_seq_x_aa,aa)
				else:
					aa="del"
					read_seq_x_aa=np.append(read_seq_x_aa,aa)
				del tripletta_nucl[:] # una volta tradotta la tripletta,cancelliamo i nucleotidi presenti e ricominciamo con i prossimi tre nucleotidi


		#-------- checkpoint!!!! ------------------
		if len(read_seq_x_aa) != len(correct_pos_ref_aa_seq_x):
			raise RuntimeError("lenght of read seq translated is different from the lenght of the pos_ref_aa array:iteration %" % (x))

		#-------------------------------------------
		 		
		read_seq_x_aa=np.array([read_seq_x_aa])		
		array_pos_ref_read_seq_x_aa=np.array([correct_pos_ref_aa_seq_x])
		array_pos_ref_read_seq_trad=np.concatenate((array_pos_ref_read_seq_x_aa,read_seq_x_aa),axis=0).T
		df_pos_ref_read_seq_trad=pd.DataFrame(array_pos_ref_read_seq_trad,columns=["pos_ref","read_aa"])
		df_pos_ref_read_seq_trad.set_index("pos_ref",inplace=True)
		lista_dataframes_pos_ref_read_trad.append(df_pos_ref_read_seq_trad)

	return lista_dataframes_pos_ref_read_trad


# funzione per costruire la tabella di prevalenza di tutte posizioni per ogni aminoacido possibile
def table_prevalence():
	# costruiamo il dataframe dove riportiamo tutte le posizioni della reference(righe)
	# e tutti i possibili aminoacidi che possono essere presenti. Riportiamo la prevalence  rispetto al totale delle reads mappate di ogni aminoacido
	Iupac="ACDEFGHIKLMNPQRSTVWYX*-" # - indica delezione,+ indica inserzione,* indica aminoacido sconosciuto o stop codon interno alla sequenza
	# * means stop codon. You can use simple methods to translate DNA to protein and they will result with internal stop codons, it's all a matter of what you defined (frame, stop at stop or continue). If you see a sequence like this it's probably not a proper protein sequence.
	Iupac=list(Iupac)
	Iupac[-1]="del"
	array_index_vs_aa_ref_of_mapping,array_index_ref_1a=reference_selection_and_trad(translation=True)
	df_index_vs_aa_ref_of_mapping=pd.DataFrame(array_index_vs_aa_ref_of_mapping,columns=["pos","ref_aa"])
	df_index_vs_aa_ref_of_mapping.set_index("pos",inplace=True)
	dataframe_all_ref_pos_all_aa=pd.DataFrame(0,index=df_index_vs_aa_ref_of_mapping.index,columns=Iupac,dtype=int)
	# creiamo il dataframe dove inserire le informazioni che si riferiscono agli aminoacidi di inserzioni
	Iupac_ins=['A+', 'C+', 'D+', 'E+', 'F+', 'G+', 'H+', 'I+', 'K+', 'L+', 'M+', 'N+', 'P+', 'Q+', 'R+', 'S+', 'T+', 'V+', 'W+', 'Y+', 'X+','*+']
	dataframe_all_ref_pos_aa_ins=pd.DataFrame(0,index=df_index_vs_aa_ref_of_mapping.index,columns=Iupac_ins,dtype=int)
	#--------- adesso riportiamo le conte nei dataframes ------------------------
	lista_dataframes_pos_ref_read_trad=translation_of_the_ref_portion_mapped()

	for i in range(len(lista_dataframes_pos_ref_read_trad)):
		array_aa_ins_based_on_pos=[]
		for aa in dataframe_all_ref_pos_all_aa.columns:
			p=0
			for pos in lista_dataframes_pos_ref_read_trad[i].index.values.astype(str):
				p+=1
				if "+" not in pos:
					try:
						if lista_dataframes_pos_ref_read_trad[i].loc[pos]["read_aa"] == aa:
							dataframe_all_ref_pos_all_aa.loc[pos][aa]+=1
 						else:
							continue
					except:
						warnings.warn("keyError ignored in dataframe_all_ref_pos_all_aa in table_prevalence module: iteration %" % (i))
						continue

		for pos in lista_dataframes_pos_ref_read_trad[i].index.values.astype(str) :
			if "+" in pos:
				# prendiamo la posizione e gli aminoacidi di inserzione e li salviamo in un array a parte
				aa_ins_pos_x=lista_dataframes_pos_ref_read_trad[i].loc[pos]["read_aa"]
				array_aa_ins_based_on_pos=np.append(array_aa_ins_based_on_pos,"%s_%s" % (pos.split("+")[0],aa_ins_pos_x))
						
		# aggiorniamo anche il dataframe che si riferisce solo alle inserzioni,aggiungendo +1 in corrispondenza della posizione dove la read ha inserzione di quell'aminoacido.
		# NOTA: NON stiamo contando QUANTI aa ha come inserzione una certa read,ma SOLO SE ha quell'aa come aminoacido di inserzione
		array_aa_ins_based_on_pos=np.unique(array_aa_ins_based_on_pos)

		for insertion in array_aa_ins_based_on_pos:
			pos_insertion=insertion.split("_")[0]
			aa_insertion=insertion.split("_")[1]
			aa_insertion=aa_insertion + "+"
			if aa_insertion in dataframe_all_ref_pos_aa_ins.columns.values:
				try:
					dataframe_all_ref_pos_aa_ins.loc[pos_insertion]["%s" % (str(aa_insertion))]+=1
				except:
					warnings.warn("keyError ignored in dataframe_all_ref_pos_aa_ins in table_prevalence module: iteration %" % (i))
					continue


	#----------- adesso sostituiamo gli indices della ref di mapping con gli indices della ref 1a ------------------------------
	# per entrambi i dataframes
	# se voglio che la table of prevalence abbia gli indici della ref of mapping e non della ref 1a cancello il seguente codice.
	array_index_ref_1a=array_index_ref_1a[0]
	dataframe_all_ref_pos_all_aa.set_index(array_index_ref_1a.astype(str),drop=True,inplace=True)
	dataframe_all_ref_pos_aa_ins.set_index(array_index_ref_1a.astype(str),drop=True,inplace=True)
	# Nota: gli indices della ref of mapping possono essere uguali a quelli della ref 1a o essere molto diversi dipende dall'allineamento
	# fatto precedentemente fra le due ref.


	#----------- calcoliamoci le prevalenze rispetto al totale delle reads mappate e analizzate ---------------------
	_,n_reads_mapped_correct_and_analyzed=sam_reads_information(output="mapped_correct")
	df_all_ref_pos_all_aa_prevalence=dataframe_all_ref_pos_all_aa.divide(n_reads_mapped_correct_and_analyzed)
	df_all_ref_pos_aa_ins_prevalence=dataframe_all_ref_pos_aa_ins.divide(n_reads_mapped_correct_and_analyzed)
	df_all_ref_pos_all_aa_prevalence=df_all_ref_pos_all_aa_prevalence.round(3)
	df_all_ref_pos_aa_ins_prevalence=df_all_ref_pos_aa_ins_prevalence.round(3)
	# se in corrispondenza di un certo aminoacido per una certa posizione abbiamo alla fine zero,significa che in quella posizione nessuna reads ha mai mappato con quel determinato aminoacido
	# nota inoltre che non è detto che la somma delle percentuali presenti in ogni posizione debba fare 100%. Magari solo una parte delle reads ha mappato in quella posizione.
	# es: 4 reads mappate. In pos 420: 1 read ha una G, 1 read ha una C, le altre 2 non hanno mappato in 420. Rientrano fra le reads che hanno mappato sulla ref,ma non hanno mappato in quel punto. 
	# Quindi su tutte le reads mappate, 25% ha in pos 420 una G,un altro 25% ha in pos 420 una C,il 50% non ha mappato in 420.
	return df_all_ref_pos_all_aa_prevalence,df_all_ref_pos_aa_ins_prevalence



rule_set_folder=os.path.abspath("/".join((os.path.abspath("/".join(sys.argv[0].split("/")[:-1]))).split("/")[:-1])) + "/rule_set"	
def compare_to_the_rule_set(export_results):
	_,_,_,n_reads_mapped_correct_and_analyzed=sam_reads_information(output="all")
	df_all_ref_pos_all_aa_prevalence,df_all_ref_pos_aa_ins_prevalence=table_prevalence()
	genotype_of_ref_of_mapping=sam_reads_information(output="genotype_ref_map")
	# apriamo il rule set e parsiamolo come dataframe
	open_input="%s/rule_set_geno2pheno.txt" % (rule_set_folder)
	file_rule_set=open(open_input,"r")
	rule_set=[x.strip().split("\t") for x in file_rule_set]
	df_rule_set=pd.DataFrame(rule_set,columns=["Drugs","Region","Rule","Subtype","Predictions","Reference"])
	# prendiamo a parte tutti farmaci presi in considerazione nel rule set
	all_drugs_rule_set=np.array([df_rule_set["Drugs"]])
	all_drugs_rule_set=np.unique(all_drugs_rule_set)
	# -------------seleziono la parte del rule set dedicata al gene che stiamo analizzando e al genotipo del nostro campione------------
	df_rule_set=df_rule_set[df_rule_set["Region"]=="%s" % (dictionary_args["gene"])]
	# ------------ seleziono la parte del rule set dove compare il genotipo del campione analizzato ------------------------------------
	# creiamo la lista dei genotipi presenti nella colonna subtype che il nostro dataframe dovrò tenere in considerazione
	# infatti non è detto che il genotipo di input sia scritto allo stesso modo di come è scritto nel file del rule set,inoltre la rule si potrebbe riferire a più genotipi in contemporanea)
	lista_genotypes_set=[] 
	for genotypes in df_rule_set["Subtype"]:
		if  genotype_of_ref_of_mapping in genotypes:
			lista_genotypes_set.append(genotypes)
	df_rule_set=df_rule_set.loc[df_rule_set["Subtype"].isin(lista_genotypes_set)]

	#------------- prendiamo le rules(pos e aa rilevante) significative del rule set ----------------------------------------------------
	# creaiamo la lista dei rule set significativi presente nel rule set divisi in base alla pos e all'aminoacido
	# infatti potrebbero esserci più rule set per lo stesso farmaco
	lista_significant_rule_set_pos=[]
	lista_significant_rule_set_aa=[]
	for rules in df_rule_set["Rule"]:
		iterator_pos=re.finditer('\d+', rules)
		iterator_aa=re.finditer('[a-zA-Z]+', rules)
		for pos in iterator_pos:
			lista_significant_rule_set_pos.append(pos.group())
		for aa in iterator_aa:
			lista_significant_rule_set_aa.append(aa.group())

	#----------- vediamo se le pos e gli aa rilevanti si trovano nella table of prevalence --------------------------------------------

	# cerchiamo le pos e gli aa rilevanti nella table of prevalence e vediamo se ci sono reads che hanno quel aa in quella pos,quindi vediamo se il valoreè diverso da 0
	lista_relevant_pos_aa_found_table_prevalence=[]
	lista_prevalence_value_pos_aa_found=[]
	p=0

	for pos in lista_significant_rule_set_pos:
		p+=1
		a=0
		for aa in lista_significant_rule_set_aa:
			a+=1
			if aa!="any" and aa!="ins":
				if df_all_ref_pos_all_aa_prevalence.loc[pos][aa] != 0.0 and a==p:
					lista_relevant_pos_aa_found_table_prevalence.append(str(pos)+str(aa))
					lista_prevalence_value_pos_aa_found.append(df_all_ref_pos_all_aa_prevalence.loc[pos][aa])
				else:
					continue
			elif aa=="any":
				Iupac="ACDEFGHIKLMNPQRSTVWYX"
				lista_Iupac=list(Iupac)
				for aa in lista_Iupac:
					if pos == lista_significant_rule_set_pos[lista_significant_rule_set_aa.index("any")]:
						if df_all_ref_pos_all_aa_prevalence.loc[pos][aa] != 0.0:
							lista_relevant_pos_aa_found_table_prevalence.append(str(pos)+"any:"+str(aa))
							lista_prevalence_value_pos_aa_found.append(df_all_ref_pos_all_aa_prevalence.loc[pos][aa])
						else:
							continue
			# questo in realtà è un codice sperimentale perchè non so come indicheranno in futuro le inserzioni nel rule set di geno2pheno
			elif aa=="ins":
				Iupac="ACDEFGHIKLMNPQRSTVWYX"
				lista_Iupac=list(Iupac)
				for aa in lista_Iupac:
					if pos == lista_significant_rule_set_pos[lista_significant_rule_set_aa.index("ins")]:
						if df_all_ref_pos_aa_ins_prevalence.loc[pos][str(aa)+"+"] != 0.0:
							lista_relevant_pos_aa_found_table_prevalence.append(str(pos)+"ins:"+str(aa))
							lista_prevalence_value_pos_aa_found.append(df_all_ref_pos_aa_ins_prevalence.loc[pos][str(aa)+"+"])
						else:
							continue

	#--------------------- dalle rules del rule set selezionato prima estraiamo solo quelle(e relative informazioni) trovate nella table of prevalence --------

	rules_selected=df_rule_set["Rule"].values # le rules che specifiche per il gene  e genotipo di interesse che abbiamo estratto dal rule set a disposizione
	# dunque estraiamo dal rule specifico per il gene e il genotipo di interesse le mutazioni rilevanti che sono state trovate nel nostro campione			
	r=-1
	array_indices_rules_extracted=np.array([])
	for rule in rules_selected:
		r+=1
		for relevant_rule_founded in lista_relevant_pos_aa_found_table_prevalence:
			if "any" not in relevant_rule_founded and "ins" not in relevant_rule_founded:
					if relevant_rule_founded in rule:
						array_indices_rules_extracted=np.append(array_indices_rules_extracted,r)
			elif "any" in relevant_rule_founded and "ins" not in relevant_rule_founded: 
					if relevant_rule_founded.split(":")[0] in rule:
						array_indices_rules_extracted=np.append(array_indices_rules_extracted,r)
			elif "any" not in relevant_rule_founded and "ins" in relevant_rule_founded:
					if relevant_rule_founded.split(":")[0] in rule:
						array_indices_rules_extracted=np.append(array_indices_rules_extracted,r)						  
	array_indices_rules_extracted=np.unique(array_indices_rules_extracted)
	# creiamo il dataframe che è una sottosezione del dataframe del rule set principale. In questo dataframe sono riportate solo le rules trovate nel campione di reads analizzato
	df_rule_set=df_rule_set.iloc[array_indices_rules_extracted]
	# creiamo il dataframe che indica le prevalence delle mutazioni trovate
	array_relevant_pos_aa_found_table_prevalence=np.array([lista_relevant_pos_aa_found_table_prevalence])
	array_prevalence_pos_aa_found=np.array([lista_prevalence_value_pos_aa_found])
	array_rules_found_vs_prevalences=np.concatenate((array_relevant_pos_aa_found_table_prevalence,array_prevalence_pos_aa_found),axis=0).T
	df_rules_found_vs_prevalences=pd.DataFrame(array_rules_found_vs_prevalences,columns=["rules_found","prevalence"])
	# ------------  selezione rules trovate in base a threshold value ------------------------------------ 
	# se è stato selezionato un threshold value,selezioniamo solo quelle rules trovate nelle reads con una prevalenza superiore alla soglia scelta
	if dictionary_args["threshold_value"]!="all_values":
		df_rules_found_vs_prevalences=df_rules_found_vs_prevalences.loc[df_rules_found_vs_prevalences["prevalence"].astype(float)>float(dictionary_args["threshold_value"])]
		rules_selected_based_on_prevalence=df_rules_found_vs_prevalences["rules_found"].values
		df_rule_set=df_rule_set.loc[df_rule_set["Rule"].isin(rules_selected_based_on_prevalence)]
	# creiamo gli arrays delle drugs resistenti e di quelle suscettibili
	array_drugs_resistant=np.array(df_rule_set["Drugs"])	
	array_drugs_susceptible=np.array(list(itertools.ifilterfalse(lambda x: x in array_drugs_resistant,all_drugs_rule_set)))
	# ---------- aggiungiamo il coverage al dataframe della table of prevalence ------------------------
	coverage_percent=df_all_ref_pos_all_aa_prevalence.iloc[:,:].sum(axis=1) # coverage in percentuale (frequenza relativa delle reads che hanno mappato in una certa posizione rispetto al totale delle reads mappate). Escluso ovviamente le inserzioni,perchè si riferiscono a posizioni aggiuntive che non esistono nella ref.
	coverage_abs=coverage_percent*n_reads_mapped_correct_and_analyzed # qua invece calcoliamo il coverage in frequenza assoluta.
	coverage_abs=coverage_abs.astype(int)
	coverage=pd.concat((coverage_percent.rename("cov_rel_freq"),coverage_abs.rename("cov_abs_freq")),axis=1)
	df_all_ref_pos_all_aa_prevalence=pd.concat((df_all_ref_pos_all_aa_prevalence,coverage),axis=1) # aggiungiamo il coverage al dataframe delle prevalence(solo quello principale,non quello dedicato alle inserzioni)
	# ---- aggiugiamo la percentuale totale di reads con inserzioni in una certa pos ------------------------------------------------------------------------------------------
	coverage_percent=df_all_ref_pos_aa_ins_prevalence.iloc[:,:].sum(axis=1)
	coverage_abs=coverage_percent*n_reads_mapped_correct_and_analyzed # qua invece calcoliamo il coverage in frequenza assoluta.
	coverage_abs=coverage_abs.astype(int)
	coverage=pd.concat((coverage_percent.rename("tot_rel_freq"),coverage_abs.rename("tot_abs_freq")),axis=1)
	# aggiungiamo la percentuale di reads mappate con inserzioni in una certa posizione al dataframe delle inserzioni(non è un vero e proprio coverage)
	df_all_ref_pos_aa_ins_prevalence=pd.concat((df_all_ref_pos_aa_ins_prevalence,coverage),axis=1)
	#--------------------- se voglio esporto fuori dallo script i risultati dell'analisi ---------------------
	if export_results==True:
		# esportiamo le tabelle di prevalenze come file csv
		df_all_ref_pos_all_aa_prevalence.to_csv("./df_all_ref_pos_all_aa_prevalence.csv")		
		df_all_ref_pos_aa_ins_prevalence.to_csv("./df_all_ref_pos_aa_ins_prevalence.csv")		
		# esportiamo i dataframe delle regole trovate con le prevalence delle mutazioni
		df_rules_found_vs_prevalences.to_csv("./df_rules_found_vs_prevalences.csv")
		df_rule_set.to_csv("./df_rule_set.csv")
		# salviamo i numpy array delle drugs resistenti e suscettibili come txt file
		with open("./drugs_resistant_and_susceptible.txt","w") as drugs_output_file:
			drugs_output_file.write("drugs_resistant:%s\ndrugs_susceptible:%s" % (array_drugs_resistant,array_drugs_susceptible))
	
	return array_drugs_resistant,array_drugs_susceptible,df_rules_found_vs_prevalences,df_rule_set

	


######################################################################################################################################
if __name__=="__main__":
	#--------------------- printiamo inanzitutto le informazioni del mapping --------------------------------
	print "---------------------------------------------------------"
	print "mapping informations:"
	n_tot_reads,n_reads_not_map,n_reads_map,n_reads_mapped_correct_and_analyzed=sam_reads_information(output="all")
	print "tot_reads_input:",n_tot_reads
	print "tot_reads_mapped_correct_and_analyzed:",n_reads_mapped_correct_and_analyzed
	print "---------------------------------------------------------"
	#-----------------  controlliamo la percentuale di reads mappate correttamente e analizzate --------------
	check_reads_information()
	#------------------- restituiamo l'output finale ---------------------------------------------------------
	df_all_ref_pos_all_aa_prevalence,df_all_ref_pos_aa_ins_prevalence=table_prevalence()
	array_drugs_resistant,array_drugs_susceptible,df_rules_found_vs_prevalences,df_rule_set=compare_to_the_rule_set(export_results=True)
	if df_rule_set.empty==True:
		raise Warning("no relevant mutations founded in the selected input reads file...")
	else:
		print "relevant mutations founded!"
		print df_rule_set

#parse_sam()
#check_reads_information()
#reference_selection_and_trad(translation=True)
#table_prevalence()
#translation_of_the_ref_portion_mapped()
#extract_read_portion_mapped()
#extract_ref_portion_mapped()
#compare_to_the_rule_set(export_results=False)
