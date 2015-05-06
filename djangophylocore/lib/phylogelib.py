#!/usr/bin/env python
#################################################
# Licence :
# Released under CeCILL license
# This program is Copyleft. You can disrtribued in
# order that you follow the CeCILL licence :
# http://www.cecill.info/
#
# Author:
# Nicolas Clairon : clairon [_at_] gmail.com
#################################################
""" Pylogenic trees manipulating fonctions """

from pyparsing import *

class NewickParser( object ):
    """
    Class wrapping a parser for building Trees from newick format strings
    """
    def __init__( self ):
        self.parser = self.create_parser()
        self.tree = []
        self.taxa_list = []
        self.source = ""

    def create_parser( self ):
        """
        Create a 'pyparsing' parser for newick format trees roughly based on the
        grammer here:
            http://evolution.genetics.washington.edu/phylip/newick_doc.html
            
        Problems:
            - Is a single leaf a valid tree?
            - Branch length on root? Doesn't make sense to me, and forces the root
              to be an edge.
        """
        # Basic tokens
        real = Combine( Word( "+-" + nums, nums ) + 
                        Optional( "." + Optional( Word( nums ) ) ) +
                        Optional( CaselessLiteral( "E" ) + Word( "+-" + nums, nums ) ) )
        lpar = Suppress( "(" ) 
        rpar = Suppress( ")" )
        colon = Suppress( ":" )
        semi = Suppress( ";" )
        quot = Suppress( "'" )
        # Labels are either unquoted or single quoted, if unquoted underscores will be replaced with spaces
        quoted_label = QuotedString( "'", None, "''" ).setParseAction( lambda s, l, t: t[0] )
        #simple_label = Word( ' 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'*+-./<=>?@[\\]^_`{|}~' )#.setParseAction( lambda s, l, t: t[0].replace( "_", " " ) )
        #VR
        simple_label = Word( ' 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'*+-./<=>?@[\\]^_`{|}~' )
        label = quoted_label | simple_label
        # Branch length is a real number (note though exponents are not in the spec!)
        branch_length = real.setParseAction( lambda s, l, t: float( t[0] ) )
        # Need to forward declare this due to circularity
        node_list = Forward()
        # A leaf node must have a label
        leaf = ( label + Optional( colon + branch_length, "" ) ) \
            .setParseAction( lambda s, l, t: t[0])
            #.setParseAction( lambda s, l, t: ( t[1] or None,  t[0]  ) )
            #.setParseAction( lambda s, l, t: Edge( t[1] or None, Tree( t[0] ) ) )
        # But a subtree doesn't but must have edges
        subtree = ( node_list + Optional( label, "" ) + Optional( colon + branch_length, "" ) )\
            .setParseAction( lambda s, l, t: [t[0]] )
            #.setParseAction( lambda s, l, t: ( t[2] or None, ( t[1] or None, t[0] ) ) )
            #.setParseAction( lambda s, l, t: Edge( t[2] or None, Tree( t[1] or None, t[0] ) ) )
        # A tree is then a nested set of leaves and subtrees
        node = leaf | subtree
        node_list << ( lpar + delimitedList( node ) + rpar ) \
            .setParseAction( lambda s, l, t: [ t.asList() ] )
        # The root cannot have a branch length
        tree = ( node_list + Optional( label, "" ) + Optional( semi ) )\
            .setParseAction( lambda s, l, t: ( t[1] or None, t[0] ) )
            #.setParseAction( lambda s, l, t: Tree( t[1] or None, t[0] ) )
        # Return the outermost element
        return tree

    def parse_string( self, s ):
        self.source = s
        self.tree = self.parser.parseString( str(s) )[0][1]
        self.taxa_list = []
        return self.tree

    def __remove_taxa( self, tree, node, remove_list ):
        if node in remove_list:
            tree.remove( node )
        else:
            if type(node) is list:
                for i in node[:]:
                    self.__remove_taxa( node, i, remove_list )

    def __remove_singleton( self, node ):
        #print ">>>node", node
        if type( node ) is list:
            for son in node:
                #print "son", son, "node", node
                if type( son ) is list and len( son ) == 1:
                    #print "before remove1", node, "<<<<<->>>>>", son
                    node[node.index( son )] = son[0]
                    #print "after remove1", node, "<<<<->>>>", son
                    self.__remove_singleton( node )
                elif type( son ) is list and len( son ) == 0:
                    #print "before remove0", node, "<<<<<->>>>>", son
                    node.remove( [] )
                    self.__remove_singleton( node )
                    #print "after remove0", node, "<<<<->>>>", son
                self.__remove_singleton( son )

    def remove_singleton( self, node ):
        """ remove singleton of the subtree rooted by node """
        if type( node ) is list:
            for son in node:
                if type( son ) is list and len( son ) <= 1:
                    self.__remove_one_singleton( node )
            for son in node:
                if type( son ) is list:
                    self.remove_singleton( son )
            self.__remove_one_singleton( node )

    def __remove_one_singleton( self, node ):
        """ remove singleton of the son of node """
        if type( node ) is list:
            for son in node:
                if type( son ) is list and len( son ) == 1:
                    node[node.index( son )] = son[0]
                    self.__remove_one_singleton( node )
                else:
                    if type( son ) is list and len(son) == 0:
                        node.remove( [] )
                        self.remove_singleton( node )

    def filter( self, remove_taxa_list ):
        tree = self.tree
        self.__remove_taxa( tree, tree, remove_taxa_list )
        self.remove_singleton( tree )
        return tree

    def __correct_tree( self, tree, node, correct_dict ):
        if type(node) is not list:
            if node in correct_dict:
                tree[tree.index( node )] = correct_dict[node]
            elif " ".join( node.split()[:2]) in correct_dict:
                old_node = node
                base_node = correct_dict[" ".join( node.split()[:2] )]
                node = base_node+" "+" ".join( node.split()[2:] )
                tree[tree.index( old_node )] = node
            elif node.split()[0] in correct_dict:
                old_node = node
                base_node = correct_dict[node.split()[0]]
                node = base_node+" "+" ".join( node.split()[1:] )
                tree[tree.index( old_node )] = node
        else:
            if type(node) is list:
                for i in node:
                    self.__correct_tree( node, i, correct_dict )

    def correct_tree( self, correct_dict ):
        """
        correct taxa name with the corresponding name in correct_dict
        {'ratus':'rattus', 'echinops':'echinops <plant>'}
        """
        tree = self.tree
        self.__correct_tree( tree, tree, correct_dict )
        return tree

     
    def get_taxa( self, tree = None ):
        if not tree:
            tree = self.tree
        self.taxa_list = str( tree ).replace( "[u'", "['").replace(", u'", ", '" )
        self.taxa_list = self.taxa_list.replace( "['", "[" ).replace( "']", "]" )
        self.taxa_list = self.taxa_list.replace( "',", "," ).replace( ", '", ", " )
        self.taxa_list = self.taxa_list.replace( '[', '' ).replace( ']', '' )
        if self.taxa_list:
            return self.taxa_list.split(', ')
        else:
            return []

    def get_struct( self ):
        """
        split the collection or tree into a list

        >>> getStruct( '((a,b),c);(((aee,bcd),(c,d)),e);' )
        ['((', 'a', ',', 'b', '),', 'c', ');(((', 'aee', ',', 'bcd', '),(', 'c', ',', 'd', ')),', 'e', ');']
        """
        struct = []
        taxa_name = []
        delim = []
        delimiter = False
        taxa = False
        for c in self.source:
            if c in '(),;':
                if taxa:
                    delimiter = True
                    taxa = False
                    struct.append( ''.join( taxa_name ) )
                    taxa_name = []
                delim.append( c )
            else:
                if not taxa:
                    delimiter = False
                    taxa = True
                    struct.append( ''.join( delim ) )
                    delim = []
                taxa_name.append( c )
        struct.append( ''.join( delim ) )
        return struct


def tidyNwk( nwk ):
    """ Strip all space and backline from newick string """
    nwk = nwk.replace( "\n", " " )
    while "  " in nwk:
        nwk = nwk.replace( "  ", " " )
    while '( ' in nwk:
        nwk = nwk.replace( '( ', '(' )
    while ' )' in nwk:
        nwk = nwk.replace( ' )', ')' )
    nwk = nwk.replace( ", (", ",(" )
    nwk = nwk.replace( " ,", "," )
    nwk = nwk.replace( ", ", "," )
    nwk = nwk.replace( ") ,", ")," )
    nwk = nwk.replace( ") )", "))" )
    nwk = nwk.replace( "( (", "((" )
    return nwk.strip()

def checkNwk( nwk ):
    """
    Check if nwk is in the correct newick format
    """
    try:
        nwk2list( nwk )
        return True
    except:
        return False
        
def removeNexusComments( nex ):
    """ remove all nexus comments """
    nwk_without_comment = ""
    for i in nex.split("["):
        if len(i.split("]")) > 1:
            nwk_without_comment += "".join( i.split("]")[1] )
        else:
            nwk_without_comment += i 
    return nwk_without_comment

def getParent(tree,node):
  """ Return the parent of the node """
  if tree == node:
    return ""
  for child in getChildren(tree):
    if child == node:
      return tree
    if node in child:
      parent = getParent(child,node)
  return parent


def getChildren(tree):
  """ Extract all first child of tree to a list """
  p = 0
  chaine = ""
  newarb = []
  for i in tree:
    if i == "(":
      p += 1
    if i == ")":
      p -= 1
    if (i == "," and p == 1) or (i == ")" and p == 0):
      newarb.append(chaine.strip())
      chaine = ""
      continue
    if i == "(" and p == 1:
      continue
    chaine += i
  return newarb

def getBrothers( tree, node ):
  """ return all node brothers in tree """
  tree = removeBootStraps( tree )
  return getChildren( getParent( tree, node ) )
 
def getNodes(tree):
  """ Return the nodes list of tree """
  listNodes = []
  if tree not in listNodes:
    listNodes.append(tree)
  for child in getChildren(tree):
    if child not in listNodes:
      listNodes.extend(getNodes(child))
  return listNodes

def getEdges(tree):
  """ Return the edges list in tree """
  if len(getChildren(tree)) == 0:
    return []
  listEdges = []
  for child in getChildren(tree):
    listEdges.append((tree,child))
    listEdges.extend(getEdges(child))
  return listEdges

def getTree(tree):
  """ Convert a newick string to a list :
     [ [nodelist] [edgeslist] ] """
  if len(getChildren(tree)) == 0:
    return [tree,[]]
  listEdges = []
  listNodes = [tree]
  for child in getChildren(tree):
    r = getTree(child)
    listNodes.extend(r[0])
    listEdges.append([tree,child])
    listEdges.extend(r[1])
  return [listNodes,listEdges]

def removeBootStraps(tree):
  """ Remove all bootstraps from tree """
  chaine = ""
  ignore = False
  for i in range(len(tree)):
    if tree[i] == ":":
      ignore = True
    if tree[i] == "," or tree[i] == ")":
      ignore = False
    if i > 0 and tree[i-1] == ')' and tree[i] not in "),":
      ignore = True
    if not ignore:
      chaine += tree[i]
  return chaine

#def getTaxa(tree):
#  """ Return the taxas list """
#  tree = removeBootStraps(tree).strip()
#  if len(getChildren(tree)) == 0:
#    return [tree]
#  l = []
#  for child in getChildren(tree):
#    l.extend(getTaxa(child))
#  return l

def getTaxa( tree ):
    return [i.strip() for i in removeBootStraps(tree.strip()).replace("(","").replace(")","").split(",")]

def getStruct( tree ):
    """
    split the collection or tree into a list

    >>> getStruct( '((a,b),c);(((aee,bcd),(c,d)),e);' )
    ['((', 'a', ',', 'b', '),', 'c', ');(((', 'aee', ',', 'bcd', '),(', 'c', ',', 'd', ')),', 'e', ');']
    """
    struct = []
    taxa_name = []
    delim = []
    delimiter = False
    taxa = False
    for c in tree:
        if c in '(),;':
            if taxa:
                delimiter = True
                taxa = False
                struct.append( ''.join( taxa_name ) )
                taxa_name = []
            delim.append( c )
        else:
            if not taxa:
                delimiter = False
                taxa = True
                struct.append( ''.join( delim ) )
                delim = []
            taxa_name.append( c )
    struct.append( ''.join( delim ) )
    return struct

def getDepth(tree):
  """ Return the depth of the tree """
  if len(getChildren(tree)) == 0:
    return 0
  depth = 0
  for child in getChildren(tree):
    depth = max(getDepth(child),depth)
  return 1 + depth

def getTriplets(tree):
  """ Extract all triplets of the tree in the list listTriplets """
  listTriplets = []
  for child in getChildren(tree):
    listTaxasChild = getTaxa(child)
    listIn = []
    for i in listTaxasChild:
      for j in listTaxasChild:
        if i!=j and ([i,j] not in listIn) and ([j,i] not in listIn):
          listIn.append([i,j])
    listExt = [i for i in getTaxa(tree) if i not in listTaxasChild]
    for i in listExt:
      for j in listIn:
        chaine = "("+i+",("+j[0]+','+j[1]+"))"
        listTriplets.append(chaine)
    listTriplets.extend(getTriplets(child))
  return listTriplets

def countTriplets(tree):
  """ Count the number of triplets in the tree """
  listNodes = getNodes(tree)
  listNodes.remove(tree)
  for i in getTaxa(tree):
    listNodes.remove(i)
  result = 0
  for child in listNodes:
    nbTaxas = len(getTaxa(child))
    result += (nbTaxas * (nbTaxas-1) * (len(getTaxa(getParent(tree,child))) - nbTaxas)) / 2
  return result

def nwk2list( nwk ):
    """
    Take a newick string and return a python list

    @nwk: string
    @return: python list
    """
    nwk = removeBootStraps( nwk )
    nwk = nwk.replace( "(", "[" ).replace( ")", "]" )
    nwk = nwk.replace( ",", "','" ).replace( "]',", "]," ).replace(",'[",",[")
    nwk = nwk.replace("[","['").replace( "['[", "[[" )
    nwk = nwk.replace( "]", "']" ).replace( "]']", "]]")
    nwk = nwk.replace("['[", "[[" )
    nwk = nwk.replace("]']", "]]" )
    return eval( nwk )

def generateTree(nbTaxas, maxChildren = 2, name = 1):
  """ Generate randomly a pylogenic tree :
      nbTaxas : number of taxa wanted
      maxChildren : number for incertitude
      name : index of the name in the alphabet """
  from random import random
  import string
  assert(nbTaxas > 0)
  assert(maxChildren > 1)
  if nbTaxas == 1:
    # intToString
    strName = ""
    while name > 0:
      strName += str(string.ascii_lowercase[(name-1) % 26])
      name = (name-1)/26
    return strName
#     return intToString(name)
  if maxChildren > nbTaxas:
    maxChildren = nbTaxas
  numChildren = int(((maxChildren - 2) + 1) * random()) + 2
  coupes = []
  coupes.append(1)
  for i in xrange(numChildren - 1):
    rand = int(((nbTaxas - 2) + 1) * random()) + 2
    while rand in coupes:
      rand = int(((nbTaxas - 2) + 1) * random()) + 2
    coupes.append(rand)
  coupes.append(nbTaxas + 1)
  coupes.sort()
  results = "("
  for i in xrange(len(coupes) - 1):
    if i != 0:
      results += ","
    results += generateTree(coupes[i + 1] - coupes[i], maxChildren, name + coupes[i] - 1)
  results += ")"
  return results


if __name__ == "__main__":
    ####################################################################################################
    """ Exemple of pylogenic trees """
    a = "((a,e),c,(d,(e,f)))"
    b = "((((a,b),c,d),(e,f)),g)"
    c = "(((a,b),c,d),e,f,((g,h),i))"
    d = "(t5:0.004647,t4:0.142159,((t6:0.142159,t1:0.047671)10:0.115,DinosauriaxDinosorus:0.545582)60:0.995353)"
    ####################################################################################################
    parser = NewickParser()
#    tree = [['mus', 'rattus'], [['homo', 'pan'], 'echinops']]
#    print parser.parse_string( '(((echinops <mammal> sp./ 2344)),petunia_x_hybrida);' )
#    print parser.get_taxa()
#    print parser.parse_string( """((batomys granti ear1822,batomys granti 458948),(archboldomys luzonensis,(chrotomys gonzalesi,(rhynchomys isarogensis,(((apomys datae 167243,apomys datae 167358),(apomys gracilostris m646,apomys gracilostris m648)),(((apomys sp. d 154816,apomys sp. d 154854),(((apomys hylocoetes 147871,apomys hylocoetes 147914,apomys hylocoetes 148149,apomys insignis 147915,apomys insignis 147924),(apomys insignis 147911,apomys insignis 148160)),(apomys sp. f 458762,apomys sp. f ear1491))),(((apomys microdon 167241,apomys microdon 167242),(apomys microdon 458907,apomys microdon 458919)),((apomys musculus 458925,apomys musculus 458913),(((apomys sp. a/c 135715,apomys sp. a/c 137024,apomys sp.  a/c 458747),apomys sp. a/c 458751),(apomys sp. b 145698,apomys sp. b 145699))))))))))""" )
#    print parser.get_taxa()
#    print parser.parse_string( '(echinops_<mammal>,petunia_x_hybrida);' )
#    print parser.get_taxa()
#    print parser.parse_string( u"(pan,petunia hybrida)" )
#    print parser.get_taxa()
#    print parser.parse_string( "(('mis france', 'rattus'), ((('homo', 'pan')), 'echinops'));" )
#    print parser.correct_tree( {"mis": "mus", "echinops":"echinops <plant>" } )
#    print parser.get_taxa()
#    print parser.parse_string( "(glis,(mus,rattus rattus));")
#    print parser.filter( ['mus', 'glis', 'rattus rattus'] )
#    print parser.get_taxa()
#    print parser.parse_string( "(mus:0.987,rattus:0.93784);")
#    print parser.get_taxa()
    #print parser.parse_string( '((((((((bla, ble))), ((bla, bli), blo),(blu, bla))))))' )
    #print "ooo", parser.filter( ['bla', 'ble'] )
    print parser.parse_string( '((bla),((bli),blo),blu)' )
    print "ooo", parser.filter( ['bla', 'ble'] )
    print parser.parse_string( '((((((((abrocomidae,chinchillidae),((agoutidae,dasyproctidae),(caviidae,hydrochaeridae)),capromyidae,dinomyidae,((echimyidae,myocastoridae),octodontidae)),erethizontidae),((bathyergidae,thryonomyidae),hystricidae)),((dipodidae,muridae),pedetidae)),(anomaluridae,(((aplodontidae,sciuridae),myoxidae),castoridae),(geomyidae,heteromyidae))),(((((((((((antilocapridae,(bovidae,(cervidae,giraffidae))),tragulidae),(((balaenidae,(balaenopteridae,eschrichtiidae)),(((((delphinidae,phocoenidae),monodontidae),platanistidae),ziphiidae),physeteridae)),hippopotamidae)),(suidae,tayassuidae)),camelidae),(equidae,(rhinocerotidae,tapiridae))),(carnivora,manidae)),((((craseonycteridae,emballonuridae),((megadermatidae,(nycteridae,rhinolophidae)),rhinopomatidae)),((((furipteridae,natalidae),thyropteridae),myzopodidae),(molossidae,vespertilionidae),(((mormoopidae,phyllostomidae),noctilionidae),mystacinidae))),pteropodidae)),((((bradypodidae,megalonychidae),dasypodidae),myrmecophagidae),((chrysochloridae,tenrecidae),((((dugongidae,trichechidae),procaviidae),elephantidae),(macroscelididae,orycteropodidae))))),(cynocephalidae,(primates,tupaiidae))),(leporidae,ochotonidae))),(erinaceidae,(solenodontidae,(soricidae,talpidae))))')
    t = parser.filter ([u'abrocomidae',
u'agoutidae',
u'ailanthus desf.',
u'anomaluridae',
u'antilocapridae',
u'aplodontidae',
u'balaenidae',
u'balaenopteridae',
u'bathyergidae',
u'black cherry',
u'bovidae',
u'bradypodidae',
u'camelidae',
u'canidae',
u'capromyidae',
u'carnivora',
u'castoridae',
u'caviidae',
u'cebidae',
u'cercopithecidae',
u'cervidae',
u'cheirogaleidae',
u'chinchillidae',
u'chrysochloridae',
u'craseonycteridae',
u'ctenodactylidae',
u'ctenomyidae',
u'cynocephalidae',
u'dasypodidae',
u'dasyproctidae',
u'daubentoniidae',
u'delphinidae',
u'dinomyidae',
u'dipodidae',
u'dugongidae',
u'echimyidae',
u'elephantidae',
u'emballonuridae',
u'equidae',
u'erethizontidae',
u'erinaceidae',
u'eschrichtiidae',
u'felidae',
u'furipteridae',
u'geomyidae',
u'giraffidae',
u'herpestidae',
u'heteromyidae',
u'hippopotamidae',
u'hominidae',
u'hyaenidae',
u'hydrochaeridae',
u'hylobatidae',
u'hystricidae',
u'lemuridae',
u'leporidae',
u'macroscelididae',
u'manidae',
u'megadermatidae',
u'megalonychidae',
u'molossidae',
u'monodontidae',
u'mormoopidae',
u'moschidae',
u'muridae',
u'mustelidae',
u'myocastoridae',
u'myrmecophagidae',
u'mystacinidae',
u'myzopodidae',
u'natalidae',
u'neobalaenidae',
u'noctilionidae',
u'nycteridae',
u'ochotonidae',
u'octodontidae',
u'odobenidae',
u'oligohymenophorans',
u'orycteropodidae',
u'otariidae',
u'pedetidae',
u'petromuridae',
u'phocidae',
u'phocoenidae',
u'photinia lindl.',
u'phyllostomidae',
u'physeteridae',
u'plasmodium knowlesi (strain nuri)',
u'platanistidae',
u'primates',
u'procaviidae',
u'procyonidae',
u'pteropodidae',
u'rhinocerotidae',
u'rhinolophidae',
u'rhinopomatidae',
u'rowan',
u'sciuridae',
u'solenodontidae',
u'soricidae',
u'suidae',
u'talpidae',
u'tapiridae',
u'tarsiidae',
u'tayassuidae',
u'tenrecidae',
u'thryonomyidae',
u'thyropteridae',
u'tragulidae',
u'trichechidae',
u'tupaiidae',
u'ursidae',
u'vespertilionidae',
u'viverridae',
u'ziphiidae'] )

    print t
