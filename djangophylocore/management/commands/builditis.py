from django.core.management.base import NoArgsCommand
from optparse import make_option

import os
import sys, csv
from django.conf import settings
import httplib, re

localDir = os.path.dirname(__file__)
absDir = os.path.join(os.getcwd(), localDir)
DUMP_PATH = os.path.join( absDir,'..','..', 'dumps' ) 
REL_ID = 0

class Command(NoArgsCommand):
    option_list = NoArgsCommand.option_list + (
        make_option('--verbose', '-v', action='store_true', dest='verbose', 
            help='Verbose operation'),
    )
    help = "Download and build the itis database"
    requires_model_validation = False

    def download_itis( self, verbose ):
        if not os.path.exists( './itisdump.tar.gz' ):
            if verbose:
                print "Downloading ITIS database on the web"
	    url_thumb = ''
	    conn = httplib.HTTPConnection("www.itis.gov")
	    conn.request("GET", "/downloads/")
	    f = conn.getresponse().read()
	    file_name = re.findall('href="[^"]*itisMS[^"]*gz',f)[0][8:]
            os.system( "curl -# http://www.itis.gov/downloads/%s > itisdump.tar.gz" % file_name)
            #os.system( "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" )
        if verbose:
            print "Extracting database... please wait"
        os.system( "tar xf itisdump.tar.gz" )
        os.system( "mv itisMS.* itisdump" )
 
    def handle_noargs(self, **options):
        global DUMP_PATH
        verbose = options.get("verbose", False)
        if verbose:
            print "loading taxonomy, please wait, it can take a while..."
        if not os.path.exists( DUMP_PATH ):
            os.system( 'mkdir %s' % DUMP_PATH )
        else:
            os.system( 'rm %s/*' % DUMP_PATH )
        self.download_itis( verbose )
        path = './itisdump'
        synonym_file_path = os.path.join( path, 'synonym_links' )
        taxonomic_units_path  = os.path.join( path, 'taxonomic_units' )
        taxon_unit_types_path  = os.path.join( path, 'taxon_unit_types' )
        vernaculars_path  = os.path.join( path, 'vernaculars' )
        kingdoms_path  = os.path.join( path, 'kingdoms' )
        # building
        # kingdoms
        d_kingdom = self.getKingdom( kingdoms_path )
        # synonyms
        syn_tax = self.getSynonym( synonym_file_path)
        # getting rank
        rank = self.getRank( taxon_unit_types_path )
        #return correct_tax, tax_name, tax_id
        #correct_tax = tax_name = tax_id = getCorrectTaxa("/Users/vranwez/Desktop/ITIS/itis_fic_utils/taxonomic_units")
        taxa_sons={}
	# collect informations on correct taxa
        correct_taxa, max_id = self.getCorrectTaxa(taxonomic_units_path, taxa_sons, syn_tax)
        root_id = max_id
	## recuperation des homonyms
        taxa_homo={}
        homonyms = {}
        for tax_id in correct_taxa:
            tax_name = correct_taxa[tax_id]["name"]
            if not tax_name in taxa_homo:
                 taxa_homo[tax_name]=[]    
            taxa_homo[tax_name].append(tax_id)
        
        
        for tax_name, tax_id in taxa_homo.items():
            if len( tax_id ) > 1:
                for id in tax_id:
                    if not tax_name in homonyms:
                        homonyms[tax_name] = []
                    # <> not close yet we may need to add itis id
                    homonyms[tax_name].append( ( id, "%s <%s %s" % (correct_taxa[id]['name'],
                      d_kingdom[correct_taxa[id]['kingdom']], rank[correct_taxa[id]['rank']] )) )
        #if some ambiguity still exist, we add itis id to the name
        taxa_homo_freq={}
        for homonym_name, scientific_names_list in homonyms.iteritems():
            for n in scientific_names_list :
                if not n[1] in taxa_homo_freq : 
                    taxa_homo_freq[n[1]]=0
                taxa_homo_freq[n[1]]= taxa_homo_freq[n[1]]+1;
                
        homonyms_tmp ={}
        for homonym_name, scientific_names_list in homonyms.iteritems():
            homonyms_tmp[homonym_name]=[]
            for n in scientific_names_list :
                if taxa_homo_freq[n[1]] > 1 : 
                    homonyms_tmp[homonym_name].append( ( n[0], "%s itis %s>" %(n[1],n[0]) ) )
                else:
                    homonyms_tmp[homonym_name].append( ( n[0], "%s>" %(n[1]) ) )
               
        homonyms = homonyms_tmp; 
          
        # avery taxon should now have an unambiguouds name
        # XXXXXXX in generate dumps the name should be find using homonyms table
        homonym_ok = {}
        for homonym_name, scientific_names_list in homonyms.iteritems():
            unambiguous_names =[n[1] for n in scientific_names_list]
            names = set(unambiguous_names)
            if len( names ) != len( unambiguous_names ): # ambiguous names even when adding kingdom to the name
                print "===========> ambiguous name should no longer exist"
                for taxa in scientific_names_list:
                    del correct_taxa[taxa[0]] # taxa[0] = id    
            else:
                homonym_ok[homonym_name] = scientific_names_list
        homonyms = homonym_ok
#		for taxa in scientific_name_list:
#                    correct_taxa[taxa[0]]["name"] = taxa[1]  # taxa[1] = new name    
        ## recuperation des taxa valid ayant un pere valid 
        reachable_taxa={}
        #202420
        taxa_homo={}
        self.compute_reachable_taxa(correct_taxa, taxa_sons, str(root_id), reachable_taxa, taxa_homo)
        print ">>>", len(reachable_taxa)
        ancestor_file = open( os.path.join( DUMP_PATH, 'parentsrelation.dmp' ), 'a')
        self.compute_write_ancestor( taxa_sons,str(root_id),[],ancestor_file)
        ancestor_file.close()
        
        # getting common names
        common_name = self.getCommonName( vernaculars_path )
        #for tax_id, tax_info in correct_tax_tree.items():
        #    parent_id = tax_info["parent_id"]
        #    print tax_id + " "+parent_id 
        # recuperation des noms de synonyms
        synonyms = {}
        taxa_names = self.getTaxaName(taxonomic_units_path)
        reachable_values={}
        for tax_id, tax_info in reachable_taxa.iteritems():
            reachable_values[tax_info["name"]] = tax_id
        
        for syn_id, ref_id in syn_tax.iteritems():
            
            syn_name = taxa_names[syn_id]
            #print "testReach:"+ref_id in reachable_taxa +"testSynNameSc"+(syn_name not in reachable_values)
            if ref_id in reachable_taxa and (syn_name not in reachable_values):
                synonyms[syn_name] = ref_id
        syn_tax = synonyms
        #generating dumps
        self.generating_dumps( max_id, reachable_taxa, rank, common_name,
          homonyms, syn_tax, taxa_sons, d_kingdom, correct_taxa, taxa_sons )
        os.system( 'iconv -f iso8859-15 -t utf-8 %s > %s' % (os.path.join(
          DUMP_PATH, 'taxonomy.dmp' ), os.path.join( DUMP_PATH, 'taxonomy.dmp_utf-8' )))
        os.system( 'mv %s %s' % ( os.path.join( DUMP_PATH,
          'taxonomy.dmp_utf-8' ), os.path.join( DUMP_PATH, 'taxonomy.dmp' ) ))

#    def clean_name(self,name):
#        cleanNameNwk = name.replace( ")", "_" ).replace( "(", "_" ).replace(",", " ").replace(":", " ").replace(";", " ")
#        cleanName = cleanNameNwk.replace("'", " ").replace("`"," ").replace('"',' ').strip()
#        return cleanName
# VR sept 09
    def clean_name(self,name):
        cleanName ="";
        allowed="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890+-*<>";
        spaceOpen = False
        for l in name :
            if l in allowed :
                cleanName = cleanName + l;
                spaceOpen = False
            else :
                if not spaceOpen :
                    spaceOpen = True
                    cleanName = cleanName+" "
        cleanName = cleanName.strip()
        #cleanName = cleanName.replace(" ", "_")
        return cleanName
        #cleanNameNwk = name.replace( ")", "_" ).replace( "(", "_" ).replace(",", " ").replace(":", " ").replace(";", " ")
        #cleanName = cleanNameNwk.replace("'", " ").replace("`"," ").replace('"',' ').strip()
              

    def generating_dumps( self, max_id, taxonomy, rank, common_name, homonyms,
      synonyms, ancestors, d_kingdom, correct_taxa, taxa_sons ):
        # rank
        open( os.path.join( DUMP_PATH, 'rank.dmp' ), 'w' ).write( '\n'.join( ['|'.join( i ) for i in rank.items()] )+'\n')
        # scientifics
        result = []
        TAXONOMY_TOC = set([])
        homonym_id = {}
        for homonym_name, scientific_names_list in homonyms.iteritems():
            for n in scientific_names_list:
                homonym_id[n[0]]=n[1];
            
        for taxa_id in taxonomy:
            if taxonomy[taxa_id]['name'] in homonyms:
               taxonomy[taxa_id]['name'] =homonym_id[taxa_id]

        for taxa_id in taxonomy:
            if not taxonomy[taxa_id]['name'] in TAXONOMY_TOC:
                result.append( "%s|%s|scientific name|%s|%s\n" % ( taxa_id,
                  taxonomy[taxa_id]['name'], taxonomy[taxa_id]['rank'], taxonomy[taxa_id]['parent_id'] ) )
                TAXONOMY_TOC.add( taxonomy[taxa_id]['name'] )
        open( os.path.join( DUMP_PATH, 'taxonomy.dmp' ), 'w' ).write( ''.join( result ) )
        # homonyms
        result_taxonomy = []
        result_rel = []
        index = 0
        for homonym_name in homonyms:
            assert homonym_name not in TAXONOMY_TOC, homonyms[homonym_name]
            max_id += 1
            #result_taxonomy.append( "%s|%s|homonym|300|\n" % ( max_id, homonym_name ) )
            result_taxonomy.append( "%s|%s|homonym|300|%s\n" % ( max_id, homonym_name,max_id ) )
            TAXONOMY_TOC.add( homonym_name )
            for taxa_id in homonyms[homonym_name]:
                index += 1
                result_rel.append( "%s|%s|%s\n" % (index, max_id, taxa_id[0] ) )
        open( os.path.join( DUMP_PATH, 'taxonomy.dmp' ), 'a' ).write( ''.join( result_taxonomy ) )
        open( os.path.join( DUMP_PATH, 'relhomonymtaxa.dmp' ), 'w' ).write( ''.join( result_rel ) )
        # synonyms
        result_taxonomy = []
        result_rel = []
        index = 0

        for synonym_name in synonyms:
            assert synonym_name not in TAXONOMY_TOC, str(synonyms[synonym_name])+" "+ synonym_name
            max_id += 1
            #result_taxonomy.append( "%s|%s|synonym|300|\n" % ( max_id, synonym_name ) )
            result_taxonomy.append( "%s|%s|synonym|300|%s\n" % ( max_id, synonym_name,max_id ) )
            TAXONOMY_TOC.add( synonym_name )
            index += 1
            result_rel.append( "%s|%s|%s\n" % (index, max_id, synonyms[synonym_name] ) )
        open( os.path.join( DUMP_PATH, 'taxonomy.dmp' ), 'a' ).write( ''.join( result_taxonomy ) )
        open( os.path.join( DUMP_PATH, 'relsynonymtaxa.dmp' ), 'w' ).write( ''.join( result_rel ) )
        # common names
        result_taxonomy = []
        result_rel = []
        index = 0
        already_done = set([])
        for name in common_name:
            nbRefOK=0
            #check that at least one corresponding scientific name is present in TOC
            for common in common_name[name]:
                if common['id'] in taxonomy:
                    nbRefOK+=1
            if name not in TAXONOMY_TOC and nbRefOK>0:
                max_id += 1
                #result_taxonomy.append( "%s|%s|common|300|\n" % ( max_id, name ) )
                result_taxonomy.append( "%s|%s|common|300|%s\n" % ( max_id, name,max_id ) )
                TAXONOMY_TOC.add( name )
                
                for common in common_name[name]:
                    index += 1
                    if not (max_id, common['id'] ) in already_done and common['id'] in taxonomy:
                        result_rel.append( "%s|%s|%s|%s\n" % (index, max_id, common['id'], common['langage'] ) )
                        already_done.add( (max_id, common['id']) )
            else:
		pass
                #print "%s is already in toc" % name
        open( os.path.join( DUMP_PATH, 'taxonomy.dmp' ), 'a' ).write( ''.join( result_taxonomy ) )
        open( os.path.join( DUMP_PATH, 'relcommontaxa.dmp' ), 'w' ).write( ''.join( result_rel ) )
 
    def getTaxaName( self, taxonomy_units_file ):
        fichier = open( taxonomy_units_file ).readlines()
        taxa_names = {}
        for line in fichier:
            ligne = line.lower().split('|')
            tax_name_base = " ".join( [ligne[2], ligne[4], ligne[6], ligne[8]]).strip()
            #tax_name = tax_name_base.replace( ")", " " ).replace( "(", " " ).replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " ")
            tax_name=self.clean_name(tax_name_base)
            tax_id = ligne[0]
            taxa_names[tax_id] = tax_name
        return taxa_names

    def getCorrectTaxa( self, taxonomic_units_file, taxa_sons, syn_tax ):
        #fichier_csv=open(taxonomic_units_file, "rb")
        #fichier = csv.reader(fichier_csv, delimiter='|')
        fichier = open( taxonomic_units_file ).readlines()
        correct_tax = {}
        #mis en 1ere ligne du fichier
        #correct_tax["0"] = {"name": "root", "parent_id": "0", "rank": ""}
        ## recuperation des taxa valid
        ## Compute root_id
        root_id = max( [int(i.split('|')[0]) for i in fichier] ) + 1
        print root_id
        for ligne in fichier:
            ligne = ligne.lower().split('|')
            tax_id = ligne[0]
            valid = ligne[10]
            tax_name_base = " ".join( [ligne[2], ligne[4], ligne[6], ligne[8]]).strip()
            #tax_name = tax_name_base.replace( ")", "_" ).replace( "(", "_" ).replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " ")
            tax_name=self.clean_name(tax_name_base)
            tax_parent_id = ligne[17]
            if tax_parent_id == "0":
                tax_parent_id = str( root_id )
            rank = ligne[21]
            kingdom = ligne[20]
            credibility_rating = ligne[12]
            if valid in ('accepted', 'valid'): # or ligne[11] in ('synonym')):
                tax_parent_id = syn_tax.get(tax_parent_id, tax_parent_id)
                # verification que pas de syno pour les taxa valid
                if tax_id in syn_tax:
                    print "valid with syno " + tax_id
                if tax_parent_id:
                    correct_tax[tax_id] = {"name": tax_name, "parent_id": tax_parent_id,
                      "rank": rank, 'kingdom': kingdom, 'credibility_rating': credibility_rating }
                    if not tax_parent_id in taxa_sons:
                        taxa_sons[tax_parent_id]=[]    
                    taxa_sons[tax_parent_id].append(tax_id)
                    #print "add  %s %s len taxa sons %i" %(tax_parent_id, tax_id, len(taxa_sons))
        #fichier_csv.close()    
        ligne = "%s||root||||||||valid||TWG standards met|unknown|unknown||1996-06-13 14:51:08|%s|||5|10|10/27/1999|" % ( root_id, root_id )
        ligne = ligne.split('|')
        correct_tax[ligne[0]] = {"name": ligne[2] , "parent_id": ligne[17],
          "rank": ligne[21], 'kingdom': ligne[20], 'credibility_rating':ligne[12]}
        return correct_tax, root_id

    def compute_reachable_taxa( self, correct_taxa, taxa_sons, taxa_id, reachable_taxa, taxa_homo):
        taxa_info = correct_taxa[taxa_id]
        reachable_taxa[taxa_id]=taxa_info
        taxa_name = taxa_info["name"];
        if not taxa_name in taxa_homo:
            taxa_homo[taxa_name]=[]    
        taxa_homo[taxa_name].append(taxa_id)
        if taxa_id in taxa_sons:
            for son_id in taxa_sons[taxa_id]:
                if son_id != taxa_id and son_id in correct_taxa: # pb de la racine 0
                    self.compute_reachable_taxa(correct_taxa,taxa_sons,son_id,reachable_taxa,taxa_homo)

    def getKingdom( self, kingdom_file ):
        fichier = open( kingdom_file ).readlines()
        kingdom = {}
        for ligne in fichier:
            ligne = ligne.lower().split('|')
            kingdom[ligne[0]] = ligne[1]
        return kingdom

    def getRank( self, rank_file ):
        #fichier_csv=open(rank_file, "rb")
        #fichier = csv.reader(fichier_csv, delimiter='|')
        fichier = open( rank_file ).readlines()
        rank = {}
        for ligne in fichier:
            ligne = ligne.lower().split('|')
            rank[ligne[1]] = ligne[2]
        #fichier_csv.close()
        return rank

    def getCommonName( self, common_name_file):
        #fichier_csv=open(common_name_file, "rb")
        #fichier = csv.reader(fichier_csv, delimiter='|')
        fichier  = open( common_name_file ).readlines()
        common_name = {}
        for ligne in fichier:
            ligne = ligne.lower().split('|')
            #print ligne[2]
            if (ligne[2] == "english" or ligne[2] == "unspecified"):
                #tax_name = ligne[1].replace( ")", " " ).replace( "(", " " ).replace(",", " ").replace(":", " ").replace(";", " ").replace("'", " ")
                tax_name = self.clean_name(ligne[1])
                if tax_name not in common_name:
                    common_name[tax_name] = []       
                common_name[tax_name].append( {"id":ligne[0],"langage":ligne[2]} )
        #fichier_csv.close()
        return common_name

  


    def getSynonym( self, syn_file):
        #fichier_csv=open(syn_file, "rb")
        #fichier = csv.reader(fichier_csv, delimiter='|')
        fichier = open( syn_file ).readlines()
        syn_tax = {}
        for ligne in fichier:
            ligne = ligne.lower().split('|')
            syn_tax[ligne[0]] = ligne[1]
      
        
        #fichier_csv.close()
        return syn_tax

    def compute_write_ancestor( self, taxa_sons, taxa_id, ancestor, ancestor_file ):
        global REL_ID
        #print "compute anc"
        index = 0
        for a in ancestor:
            REL_ID += 1
            ancestor_file.write("%s|%s|%s|%s\n" % ( REL_ID, taxa_id, a, index ))
            index += 1
        new_anc = ancestor[:]
        new_anc.insert( 0, taxa_id )
        if taxa_id in taxa_sons:
            for son_id in taxa_sons[taxa_id]:
                if son_id != taxa_id: # pb de la racine 0
                    self.compute_write_ancestor( taxa_sons, son_id, new_anc, ancestor_file )
            
    


