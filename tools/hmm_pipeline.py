# Performs HMMER/SAM model mixtures
# Edward Liaw 5/7/11

import sys
from subprocess import Popen

# Program paths
HMMBUILD  = $(BIN)/hmmbuild
HMMSCORE  = $(BIN)/hmmsearch
SAMBUILD  = $(BIN)/modelfromalign
SAMSCORE  = $(BIN)/hmmscore
HMMCONV   = $(BIN)/hmmconvert
SAMCONV   = $(BIN)/convert
LOGO      = $(BIN)/makelogo
PDF       = ps2pdf12

# Options
options = {
	"hmmer": {
		"build": {"--informat": "afa"},
		"search": {"--incdomE": 0.1}
	},
	"sam": {
		"build": {},
		"search": {"-select_mdalign":4, "-select_seq": 4, "-select_align": 4, "-select_score": 8, "-sort": 4, "-Emax": 10, "-mdEmax": 1}
	}
}
class HMMPipeLine(object):
	def __init__(self, reference_data, search_data, path=""):
		self.reference_data = reference_data
		self.search_data = search_data

	def hmmer_train():
		hmmer = HMM("hmmer", self.reference_data, path=path)
		sam   = HMM("sam", self.reference_data, path=path)
		sam_hmmer = sam.convert()
		hmmer_sam = hmmer.convert()

		hmmer.search(self.reference_data)
		sam_hmmer.search(self.reference_data)
		sam.search(self.reference_data)
		hmmer_sam.search(self.reference_data)

class HMM(object):
	def __init__(self, name, reference_data, model_file=None, path=""):
		"""Train an HMM"""
		self.name = name
		self.reference_data = reference_data
		self.model_file = model_file

		#Train the HMM if no model file is given
		if self.model_file is None:
			if name.endswith("hmmer"):
				#$(HMMBUILD) $(HMMFORMAT) $@ $(DATA)/$*$(FORMAT)
				self.model_file = "{}_{}.hmm".format(self.name, reference_data)
				options = [str(item) for opt in options["hmmer"]["build"].items() for item in opt]
				process = Popen([os.path.join(self.path, "hmmbuild"), *options, self.model_file, self.reference_data])
			elif name.endswith("sam"):
				#$(SAMBUILD) $*_sam -alignfile $(DATA)/$*$(FORMAT)
				self.model_file = "{}_{}.mod".format(self.name, reference_data)
				options = " ".join([str(item) for opt in options["sam"]["build"].items() for item in opt])
				process = Popen([os.path.join(self.path, "modelfromalign"), self.name, "-alignfile", self.reference_data, *options])
			else:
				raise RuntimeErrot("Invalid hmm name")

		process.communicate()

	def search(search_data):
		""""""
		if self.name.endswith("hmmer"):
			#$(HMMSCORE) $(HMMOPTS) -o $@.out $^ $(APIDB)
			self.results_file = "{}_{}.out"
			options = [str(item) for opt in options["hmmer"]["search"].items() for item in opt]
			process = Popen([os.path.join(self.path, "hmmsearch"), *options, "-o", self.results_file, self.model_file, search_data])
		elif self.name.endswith("sam"):
			#$(SAMSCORE) $@ -i $^ -db $(APIDB) $(SAMOPTS)
			self.results_file = "{}_{}.out"
			options = [str(item) for opt in options["hmmer"]["search"].items() for item in opt]
			process = Popen([os.path.join(self.path, "hmmsearch"), self.name, "-i", self.model_file, "-db", search_data, *options])

		process.communicate()

	def convert(self):
		if self.name.endswith("hmmer"):
			#Convert hmmer to sam
			new_name = "{}_hmmer".format(self.name)
			new_model = "{}_{}.mod".format(self.reference_data, new_name)
			hmmer2_file = "{}.hmm".format(new_name)
			with open(hmmer2_file) as hmmer2:
				process = Popen([os.path.join(self.path, "hmmconvert"), "-2", self.model_file], stdout=hmmer2)
				process.wait()
			process = Popen([os.path.join(self.path, "convert"), hmmer2_file])
			process.wait()
			os.remove(hmmer2_file)
			#SAM convert renames file to have .con.asc.mod extension
			os.rename("{}.con.asc.mod".format(new_name), new_model)
		elif self.name.endswith("sam"):
			#Convert sam to hmmer
			new_name = "{}_sam".format(self.name)
			new_model = "{}_{}.hmm".format(self.reference_data, new_name)
			# Conversion script requires .asc.mod extension
			from shutil import copy
			sam_rename = "{}_{}.asc.mod".format(self.reference_data, new_name)
			copy(self.model_file, sam_rename)
			process = Popen([os.path.join(self.path, "convert"), sam_rename])
			process.wait()
			os.rename(sam_rename)
			os.rename("{}.con.asc.hmm".format(new_name), new_model)

		return HMM(new_name, self.reference_data, model_file=new_model, path=self.path)





"""
target:
	echo ${MAKE} $$(TARGET)\_dir
	${MAKE} $$(TARGET)\_dir

%_dir:
	echo "include ${CURDIR}/Makefile" > $*/Makefile	
	basename $* | xargs -0 -I name mkdir name
	cd $*;	${MAKE} $*_all
	
%_hmmer_hmmer:	%_hmmer.hmm
	$(HMMSCORE) $(HMMOPTS) -o $@.out $^ $(APIDB)
	
%_hmmer_sam_nocal:	%_hmmer.hmm
	# Convert to HMMER 2.0 first
	$(HMMCONV) -2 $^ > $@.hmm
	$(SAMCONV) $@.hmm
	rm $@.hmm
	# Rename converted model
	mv $@.con.asc.mod $@.mod
	
	$(SAMSCORE) $@ -i $@.mod -db $(APIDB) $(SAMOPTS)
	
%_hmmer_sam:	%_hmmer.hmm
	# Convert to HMMER 2.0 first
	$(HMMCONV) -2 $^ > $@.hmm
	$(SAMCONV) $@.hmm
	rm $@.hmm
	# Rename converted model
	mv $@.con.asc.mod $@.mod
	${MAKE} $@.mlib
	
	$(SAMSCORE) $@ -i $@.mlib -db $(APIDB) $(SAMOPTS)
	
	
# Calibrated SAM model
%.mlib:	%.mod
	$(SAMSCORE) $* -modelfile $^ -calibrate 1
"""
