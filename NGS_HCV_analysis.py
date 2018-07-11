#coding:utf-8
# questi sono moduli bult in quindi sono presenti nella versione base di python
import os
import sys
import argparse
import pickle
import warnings
import pkgutil
from datetime import datetime
parser = argparse.ArgumentParser(description="HCV_software_NGS_sequences_analisis help table")
parser.add_argument("--working_directory", '-wd',default="./",help="change working directory path,default is ./,working directory is the directory from which the software is executed") # argomento per stabilire la working directory rispetto a cui gli scripts controllati da questo script avranno come riferimento. 
parser.add_argument("--reads_folder", "-R",help="reads folder path",required=True) 
parser.add_argument("--reads_processed",'-r',default=None, help="a single reads_file to be processed for mapping")
parser.add_argument("--reads_forward",'-rf', default=None,help="reads forward in a pair of reads")
parser.add_argument("--reads_reversed",'-rrv',default=None,help="reads reversed in a pair of reads")
parser.add_argument("--batching","-btch",default=None,help="analize only a random batch of n elements from the reads file:usefull for testing;default=None")
parser.add_argument("--gene", '-gn', choices=["NS3","NS5B","NS5A"],help="the gene represented by the reads",required=True)
parser.add_argument("--genotype",'-gt', choices=["1a","1b","2","3","4","unknown"], default="unknown", help="if known: insert the genotype next to the correct tag;if unknown: insert unknown next to the correct tag, blast will find the genotype of the sample; default=unknown")
parser.add_argument("--threshold_value", "-tv",default ="all_values",help="threshold value to select prevalence's values of the rules founded greater than the threshold value;default:all_values")
parser.add_argument("--paired_end_mapping","-p",action='store_true',help="specify this argument if you want perform a paired end mapping with the reads indicated by the -rf and -rrv arguments;default=False")
parser.add_argument("--processors","-prc",default=1,help="specifiy the number of processors used by blast,mothur and bwa;default=1")
parser.add_argument("--user_ref","-usrf", default=None,help="mapping against a reference chosen by the user;default=None")
args = parser.parse_args()

# -------correggiamo l'input di sintassi dello user se possibile ----------------
# prima di mettere gli argomenti nel dizionario ed esportarli correggiamo gli input che possiamo correggere.
if args.reads_folder.endswith("/"):
		args.reads_folder=args.reads_folder[:-1]
if args.reads_processed!=None:
	# nel caso in cui termina con un \
	if args.reads_processed.endswith("/"):
		args.reads_processed=args.reads_processed[:-1]
	# nel caso in cui rimangono ancora \ significa che originariamente terminava con più di una \
	if args.reads_processed.count("/",-3)>0:
		raise Exception("InputError: invalid // at the end of the reads_processed argument")
	# nel caso in cui ci sono più di una \ all'inizio
	if args.reads_processed.count("/",0,3)>1:
		raise Exception("InputError: invalid // at the start of the reads_processed argument")
	 
if args.reads_forward!=None:
	if args.reads_forward.endswith("/"):
		args.reads_forward=args.reads_forward[:-1]
	if args.reads_forward.count("/",-3)>0:
		raise Exception("InputError: invalid // at the end of the reads_forward argument ")
	if args.reads_forward.count("/",0,3)>1:
		raise Exception("InputError: invalid // at start of the reads_forward argument ")

if args.reads_reversed!=None:
	if args.reads_reversed.endswith("/"):
		args.reads_reversed=args.reads_reverse[:-1]
	if args.reads_reversed.count("/",-3)>0:
		raise Exception("InputError: invalid // at the end of the reads_reversed argument")
	if args.reads_reversed.count("/",0,3)>1:
		raise Exception("InputError: invalid // at the start of the reads_reversed argument")

if args.user_ref!=None:
	if args.user_ref.endswith("/"):
		args.user_ref=args.user_ref[:-1]
	if args.user_ref.count("/",-3)>0:
		raise Exception("InputError: invalid // at end of the user_ref argument")
	if args.user_ref.count("/",0,3)>1:
		raise Exception("InputError: invalid // at start of the user_ref argument")

#---------- esportiamo il dizionario degli argomenti corretto --------------------------------	
dictionary_args=vars(args)


#----------------------------------- controllo libraries ---------------------------------------------------------------------- 
# controlliamo che i pacchetti di python necessari a far partire il software siano installati nel pc dello user
lista_moduli_installati=[tup[1] for tup in pkgutil.iter_modules()]
if "pandas" not in lista_moduli_installati:
	raise Exception("pandas library is not installed in your python interpreter")
if "numpy" not in lista_moduli_installati:
	raise Exception("numpy library is not installed in your python interpreter")
if "Bio" not in lista_moduli_installati:
	raise Exception("Bio library is not installed in your python interpreter") 

# controlliamo che i pacchetti di bash necessari a far partire il software siano installati nel pc dello user
os_input_check_bwa_library="dpkg-query -W bwa 2>&1 | tee log_file.txt" # niente -a perchè dobbiamo sovrascrivere il log_precedente 
os.system(os_input_check_bwa_library)
os_input_check_blast_library="dpkg-query -W ncbi-blast+ 2>&1 | tee -a log_file.txt" # -a perchè aggiunge questo controllo di libreria al precedente controllo di libreria
os.system(os_input_check_blast_library)
os_input_check_mothur_library="dpkg-query -W mothur 2>&1 | tee -a log_file.txt" # -a perchè aggiunge questo controllo di libreria al precedente controllo di libreria
os.system(os_input_check_mothur_library)
with open("./log_file.txt") as log_file:
	for line in log_file:
		if "dpkg-query:" in str(line) and "bwa" in str(line):
			raise Exception("bwa bash library absent from your system,impossible to proceed")
		if "dpkg-query:" in str(line) and "ncbi-blast+" in str(line):
			raise Exception("ncbi-blast+ bash library absent from your system,impossible to proceed")
 		if "dpkg-query:" in str(line) and "mothur" in str(line):
			raise Exception("mothur bash library absent from your system,impossible to proceed")

	
#---------------------------------------------- Exception handling in caso di errato input dello user ----------------------------------------------- 
assert args.reads_folder != None, "InputError: -R argument not specified,please insert a correct reads_folder_path"
if args.reads_processed==None and args.reads_forward==None and args.reads_reversed==None:
	raise Exception("InputError: -r argument or -rf and -rrv arguments not specified, please select the reads_file to be processed or the two reads_file(forward and reverse) to be processed included inside the reads_folder")
if args.reads_processed!=None and args.reads_forward!=None and args.reads_reversed!=None:
	raise Exception("InputError:r and rf or rrv argument cannot be specified together")
if args.reads_processed!=None and (args.reads_forward!=None or args.reads_reversed!=None):
	raise Exception("InputError:r and rf or rrv argument cannot be specified together")
if args.paired_end_mapping==True and args.reads_processed!=None:
	raise Exception("InputError:-p argument and -r argument cannot be specified together")

assert args.gene!=None,"InputError: -gn argument not defined. You must specify the gene represented by the reads. The mapping of the reads will be performed on the associated gene reference"

if args.genotype=="unknown":
	warnings.warn("if you don't specify a genotype argument,the default is unknown. A blast is performed to find the virus genotype associated to the reads",Warning)

if args.threshold_value!="all_values":
	try:
		float(args.threshold_value)
	except ValueError:
		raise ValueError("InputError:the threshold chosen is not a float")

if args.reads_processed!=None:
	assert args.reads_processed.endswith(".fastq") or args.reads_processed.endswith(".gz"),"InputError: the input reads files is not in fastq format"
elif args.reads_processed==None:
	assert args.reads_forward.endswith(".fastq") or args.reads_reversed.endswith(".fastq") or args.reads_forward.endswith(".gz") or args.reads_reversed.endswith(".gz"),"InputError: the input reads files is not in .fastq or .gz format"

if args.batching != None:
	try: 
		int(args.batching)
	except:
		raise ValueError("the value of batching must be an integer")
if args.user_ref != None:
	assert args.user_ref.endswith(".fasta"), "the reference provided by the user must be in fasta format"


#-------------------------------- specifico scripts_folder_path che si trova dentro la cartella del software ---------------------------------------------------
scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) + "/scripts"

#------------------------------- exporto argomenti in un dizionario nella cartella degli scripts ---------------------------------------
pickle.dump(dictionary_args, open("%s/dictionary_args" % (scripts_folder), "w" ) )

#-------------------------------- working directory su cui verrà eseguito lo script -------------------------------------------------------
if args.working_directory != "./":
	os.chdir(args.working_directory) 
# os.chdir serve per cambiare la working directory del subprocesso,in questo caso il subprocesso è l'esecuzione di questo script,
# ma non cambia la directory del processo principale ovvero il processo che esegue il terminale. 
# Quando finisce l'esecuzione dello script e quindi degli scripts ad esso connessi si torna alla directory precedente 
	
print "current_process_working_directory_path:%s" % (os.getcwd())


#---------------------------------- fase di estrazione e checking ---------------------------------------------------------


# per l'esecuzione del primo script non metto -a, così sovrascrive il contenuto del precedente file di log con l'inizio del nuovo file log.
# Nota: il contenuto del log file riguardante il controllo delle librerie viene cancellato perchè non ci interessa più a questo punto del software
os_input="python %s/samples_extraction_and_checking.py | tee log_file.txt" % (scripts_folder)   
os.system(os_input)


#---------------------------------- fase di selezione reads e/o merging----------------------------------------------------
# selezioniamo le reads che subiranno il mapping,che derivino dal merging oppure no

if args.reads_processed==None:

	os_input="python %s/reads_selection_and_merging.py 2>&1 | tee -a log_file.txt" % (scripts_folder)  
	os.system(os_input)

elif args.reads_processed!=None:

	os_input="python %s/reads_selection_and_merging.py 2>&1 | tee -a log_file.txt" % (scripts_folder)
	os.system(os_input)

with open("./log_file.txt") as log_file:
	for line in log_file:
		if "Error" in line:
			raise RuntimeError("Error found during scripts execution")


#------------------------------- fase di blast delle reads per trovare genotipo ----------------------------------------------------------

if args.genotype=="unknown":
	#blasto le reads contro il database
	os_input="python %s/blast_and_genotype_reads.py 2>&1 | tee -a log_file.txt" % (scripts_folder) 
	os.system(os_input)

with open("./log_file.txt") as log_file:
	for line in log_file:
		if "Error" in line:
			raise RuntimeError("Error found during scripts execution")

#--------------------------------------- fase di mapping contro reference --------------------------------------------------

#Mapping delle reads sulla reference

os_input="python %s/bwa_mapping.py 2>&1 | tee -a log_file.txt" % (scripts_folder)
os.system(os_input)

with open("./log_file.txt") as log_file:
	for line in log_file:
		if "Error" in line:
			raise RuntimeError("Error found during scripts execution")

#--------------------------------------- fase di processamento del SAM: output finale --------------------------------------
#processamento SAM

os_input="python %s/processamento_SAM.py 2>&1 | tee -a log_file.txt" % (scripts_folder)
os.system(os_input)

with open("./log_file.txt") as log_file:
	for line in log_file:
		if "Error" in line:
			raise RuntimeError("Error found during scripts execution")

#------------------------------------- eliminiamo i file inutili ---------------------------
# eliminiamo le versioni parsate dei file originali,il file d merging se presente,e il blastout. Lasciamo il sam grezzo fra gli outputs.


os_input="rm -f *parsed.fasta; rm -f *trim.contigs*;rm -f blastout*;rm -f needle.txt"		
os.system(os_input)



