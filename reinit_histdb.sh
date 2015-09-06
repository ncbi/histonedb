#!/bin/bash
#This script runs all the commands to initialize hisoneDB from scratch.
echo "Reinitializing HistDB from scratch"
echo "will use a dummy file for the database, comment this out in reinit_histdb.sh if you want a full db"

echo ">gi|68448473|ref|NP_001020347.1| histone [Danio rerio]
MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQR
LVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA

>gi|21426823|ref|NP_085112.1| histone H1.1 [Mus musculus]
MSETAPVAQAASTATEKPAAAKKTKKPAKAAAPRKKPAGPSVSELIVQAVSSSKERSGVSLAALKKSLAA
AGYDVEKNNSRIKLGLKSLVNKGTLVQTKGTGAAGSFKLNKKAESKAITTKVSVKAKASGAAKKPKKTAG
AAAKKTVKTPKKPKKPAVSKKTSKSPKKPKVVKAKKVAKSPAKAKAVKPKASKAKVTKPKTPAKPKKAAP
KKK

>gi|686630953|ref|XP_009307922.1| histone H2A [Trypanosoma grayi]
MSLTGEDPLQQNPMMGPGSATADQTSIVSGGKHGGKATAARGKGKGKGKGKRGGKTGGKAGKRDKMSRAA
RADLNFPVGRIHSRLKDGLNRKQRCGASAAIYCAALLEYLTSEVIELAGAAAKTQKTERIKPRHLLLAIR
GDEELNQIVNATIARGGVVPFVHKSLEKKIIKKSKRAS

>gi|673921668|ref|NP_001288376.1| histone H2A [Zea mays]
MDSTGTGAGGKGKKGAAGRKVGGPRKKSVSRSVKAGLQFPVGRIGRYLKKGRYAQRVGTGAPVYLAAVLE
YLAAEVLELAGNAARDNKKTRIIPRHVLLAIRNDEELGKLLGGVTIAHGGVLPNINPVLLPKKTAEKASS
VGSKEAKSPKKAAKSPKKA


" > nr
find static/browse/seeds -name "*.gff" | xargs rm
find static/browse/seeds -name "*.pdf" | xargs rm

python manage.py flush << EOF
yes
EOF

#python manage.py migrate
#python manage.py sqlclear browse | python manage.py dbshell
#python manage.py sqlclear djangophylocore | python manage.py dbshell
python manage.py migrate


python manage.py buildncbi
python manage.py loadtaxonomy
python manage.py buildtaxonomytoc
python manage.py buildvariants -f 
python manage.py buildtrees -f
python manage.py buildsunburst -f
python manage.py buildblastdb -f
python manage.py buildvariantinfo -f
python manage.py buildseedinfo -f
