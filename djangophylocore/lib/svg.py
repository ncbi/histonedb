#!/usr/bin/env python
"""\
SVG.py - Construct/display SVG scenes.

The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it.

This program uses ImageMagick to display the SVG files. ImageMagick also 
does a remarkable job of converting SVG files into other formats.
"""

import os
display_prog = 'display' # Command to execute to display images.
      
class Scene:
    def __init__(self,name="svg",height=400,width=400):
        self.name = name
        self.items = []
        self.height = height
        self.width = width
        return

    def add(self,item): self.items.append(item)

    def strarray(self):
        var = ["<?xml version=\"1.0\"?>\n",
               "<svg height=\"%d\" width=\"%d\" >\n" % (self.height,self.width),
               " <g style=\"fill-opacity:1.0; stroke:black;\n",
               "  stroke-width:1;\">\n"]
        for item in self.items: var += item.strarray()            
        var += [" </g>\n</svg>\n"]
        return var


    def ensure_dir(self,f):
        d = os.path.dirname(f)
        if not os.path.exists(d):
            os.makedirs(d)
       
    def write_svg(self,filename=None, convert="convert", format="png", path=''):
        if filename:
            self.svgname = filename
        else:
            self.svgname = self.name + ".svg"
        fullFileName = os.path.join(path, self.svgname)
        self.ensure_dir(fullFileName)
        file = open (fullFileName,'w')
        file.writelines(self.strarray())
        file.close()
        if convert:
            os.system("%s %s %s" % (convert, os.path.join(path, self.svgname),os.path.join(path, "%s.%s" % (self.name, format) )))
            os.system("rm %s" % os.path.join(path, self.svgname))

    def display(self,prog=display_prog):
        os.system("%s %s" % (prog,self.svgname))
        return        
        

class Line:
    def __init__(self,start,end):
        self.start = start #xy tuple
        self.end = end     #xy tuple
        return

    def strarray(self):
        return ["  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" />\n" %\
                (self.start[0],self.start[1],self.end[0],self.end[1])]


class Circle:
    def __init__(self,center,radius,color):
        self.center = center #xy tuple
        self.radius = radius #xy tuple
        self.color = color   #rgb tuple in range(0,256)
        return

    def strarray(self):
        return ["  <circle cx=\"%d\" cy=\"%d\" r=\"%d\"\n" %\
                (self.center[0],self.center[1],self.radius),
                "    style=\"fill:%s;\"  />\n" % colorstr(self.color)]

class Rectangle:
    def __init__(self,origin,height,width,color):
        self.origin = origin
        self.height = height
        self.width = width
        self.color = color
        return

    def strarray(self):
        return ["  <rect x=\"%d\" y=\"%d\" height=\"%d\"\n" %\
                (self.origin[0],self.origin[1],self.height),
                "    width=\"%d\" style=\"fill:%s;\" />\n" %\
                (self.width,colorstr(self.color))]

class Text:
    def __init__(self,origin,text,size=24):
        self.origin = origin
        self.text = text
        self.size = size
        return

    def strarray(self):
        return ["  <text x=\"%d\" y=\"%d\" font-size=\"%d\">\n" %\
                (self.origin[0],self.origin[1],self.size),
                "   %s\n" % self.text,
                "  </text>\n"]
        
    
def colorstr(rgb): return "#%x%x%x" % (rgb[0]/16,rgb[1]/16,rgb[2]/16)


x = [1,0]*100
y = [0,1]*100
z = [0,0]*100
matrix = [
    x, y, z
]

def test():
    nb_tree = 3
    nb_taxa = 200
    scene = Scene('test', (nb_tree+1)*10, (nb_taxa+1)*10)
    pix = 10 
    j = 0
    for line in matrix:
        j += pix
        i = 0
        for col in line:
            if col == 1:
                scene.add(Rectangle((i,j),pix,pix,(255,255,255)))
            else:
                scene.add(Rectangle((i,j),pix,pix,(0,0,0)))
            i += pix 
    scene.write_svg()
    #scene.display()
    return

matrix = {1: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    2759: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    6072: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    7711: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    7742: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    7776: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    8287: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9347: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9443: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9526: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9596: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9604: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    9605: {8077: 1, 8078: 1, 8079: 0, 8080: 1},
    9606: {8077: 0, 8078: 0, 8079: 0, 8080: 1},
    9845: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    9895: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    9903: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    9989: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    10066: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    10088: {8077: 0, 8078: 0, 8079: 1, 8080: 1},
    10114: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    27592: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    32523: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    32524: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    32525: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33154: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33208: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33213: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33316: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33511: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    33553: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    35500: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    39107: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    40674: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    89593: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    91561: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    117570: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    117571: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    131567: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    207598: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    314145: {8077: 0, 8078: 1, 8079: 0, 8080: 0},
    314146: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    314147: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    314293: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    314295: {8077: 1, 8078: 1, 8079: 1, 8080: 1},
    337687: {8077: 1, 8078: 0, 8079: 1, 8080: 1},
    376913: {8077: 1, 8078: 1, 8079: 1, 8080: 1}}

def get_matrix():
    global matrix
    nb_trees = 4
    nb_taxa = 47
    scene = Scene('test', (nb_taxa+1)*10, (nb_trees+1)*10)
    pix = 10
    j = 0
    for taxa,tmp in matrix.items():
        j += pix
        i = 0
        for tree,val in tmp.items():
            if val == 1:
                scene.add(Rectangle((i,j),pix,pix,(255,255,255)))
            else:
                scene.add(Rectangle((i,j),pix,pix,(0,0,0)))
            i += pix
    scene.write_svg()
    #scene.display()
    return


if __name__ == '__main__': get_matrix()
