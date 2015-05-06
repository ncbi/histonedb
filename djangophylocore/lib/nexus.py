#-*- coding: utf-8 -*-

import phylogelib

def remove_nexus_comments( nex ):
    """ remove all nexus comments lines """
    nwk_without_comment = ""
    for i in nex.split("["):
        if len(i.split("]")) > 1:
            nwk_without_comment += "".join( i.split("]")[1] )
        else:
            nwk_without_comment += i 
    return nwk_without_comment


class Nexus( object ):

    def __init__( self, nexus_string ):
        # Nexus collection
        if not nexus_string.lower().strip()[:6] == "#nexus" or "end;" not in nexus_string.lower():
            raise ValueError, "Wrong nexus string"
        self.__nexus_ori = nexus_string
#        self.bad_trees_names = self.get_bad_trees_names()
#        if bad_trees_names:
#            raise ValueError,\
#              "The following trees are not in newick format: %s" % str(bad_trees_names)
        self.__collection = {}
        # XXX see if nexus can be uppercase
        if "translate" in nexus_string.lower():
            raw_collection = self.translate( nexus_string )
        else:
            raw_collection = nexus_string
        raw_collection = raw_collection.split( 
          'begin trees;')[1].split( 'end;')[0]
        for raw_tree in raw_collection.split( ';' ):
            if raw_tree.strip() and raw_tree.strip() != 'end':
                raw_tree = raw_tree.strip()
                tree_name, tree, rooted = self.__split_nexus_line( raw_tree )
                # tree = phylogelib.tidyNwk( tree )
                #self.__collection[tree_name] = newick.parse( tree )
                self.__collection[tree_name] = ( tree, rooted )
        
    def get_collection( self ):
        return self.__collection

    collection = property( get_collection )

    def get_bad_trees_names( self ):
        """
        Check if all trees are *really* newick trees
        """
        bad_trees_names = []
        for raw_tree in self.__nexus_ori.split('begin trees;')[1].split('end;')[0].split(';')[1:]:
            if raw_tree.strip():
                raw_tree = remove_nexus_comments( raw_tree )
                raw_tree = raw_tree.strip()
                tree_name, tree, rooted = self.__split_nexus_line( raw_tree )
                #tree = phylogelib.removeBootStraps( tree )
                if not phylogelib.checkNwk( tree.strip() ):
                    bad_trees_names.append( tree_name )
        return bad_trees_names

    def __split_nexus_line( self, nexus_tree_line ):
        """
        return a tuple with the name of the tree,  the tree in newick format
        and a boolean wich indicate if this is a rooted tree or not
        """
        tree_name, tree = nexus_tree_line.split('=')
        tree_name = tree_name.split()[1]
        tree_split = tree.split()
        if "[" in tree_split[0]:
            rooted = tree_split[0]
            tree = " ".join( tree_split[1:] )
        else:
            rooted = ""
            tree = " ".join( tree_split )
        if tree_name[0] in ["'", '"']:
            tree_name = tree_name[1:-1]
        tree = tree.strip()
        if rooted.strip() == '[&R]':
            rooted = True
        else:
            rooted = False
        return tree_name, tree, rooted
 
    def translate( self, nexus_string ):
        # Support of the nexus translate
        trans = []
        d_translation = {}
        translation = nexus_string.split("translate")[1].split(";")[0].split()
        collection = nexus_string.split("translate")[1].split(";")[1:]
        for i in range(0,len(translation),2):
            d_translation[translation[i]] = translation[i+1].replace(",","")
            trans.append( 
              (translation[i], translation[i+1].replace(',',''))
            )
        trans.reverse()
        new_collection = []
        for tree in collection:
            tree = tree.strip()
            if tree:
                if tree == 'end': break
                for indice, translation in trans:
                    if translation[0] in ["'",'"']:
                        translation = translation[1:-1]
                    tree = tree.split('=')[0]+'= '+tree.split('=')[1].replace(
                      indice, translation ).strip()
                new_collection.append( '\t'+tree.strip()+';' )
        nexus_string = "#nexus\n\nbegin trees;\n\n"+ \
          "\n".join(new_collection)+"\n\nend;"
        return nexus_string

    def get_nexus_comments( self ):
        """ remove all nexus comments """
        nwk_comment = ""
        for i in self.nexus_ori.split("["):
            if len(i.split("]")) > 1:
                nwk_comment += "".join( i.split("]")[0] )
            else:
                nwk_comment += i 
        return nwk_comment



