#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Import necessary module
import sys, os

#-----------------------------ANALIZED FILE ON FIRST PARAMETER

reads={}   #Create dictionary that containt results we use
compt1=0   #Define the counter for number total of reads
#WARNING : Si pas de fichier entré en argument ou qu'il n'existe pas /est pas trouvé (en paramètre l'utilisateur doit mettre un CHEMIN!!!!??? ?
#Check user put a file and not a directory to attemp to analyse several file on the same time
if len(sys.argv)>1 :
    if os.path.isdir(sys.argv[1]) :
        print("Ce script ne peux pas traité un dossier mais un fichier.")
    else :                                                                 #Check it's a SAM file
        entree=sys.argv[1]
        fichierEntré=entree.split("/")[-1]       #Isolate the file
        extFichier=fichierEntré.split(".")[-1]   #Isolate extension
        if extFichier!="sam" and extFichier!="SAM" : 
            print("Ce script est destiner à des fichier de type SAM.")
        else :                                                             #Check file it's not empty
            if os.path.getsize(sys.argv[1])==0 :     # PISTE AM2LIORATION : autre commande pour vérifier (empty()) ou afficher le premier read à l'utilisateur et lui demander si ca comvient
                print("Le fichier à analyser est vide. Le script ne peut pas s'exécuter.")
            else :
                file=open(sys.argv[1],"r")
                lines=file.readlines() #All file's text line per line in variable "lines" with readlines() function
                for line in lines :
                    temp=line.split("\t") #Variable temp change between each for
                                       #Délimiter for each information is the tabulation in SAM file
                    if temp[0][0]!='@' : #If the ligne don't begin with @ (@ for informative lines)  #[0][0]=première colonne premier caractère
                        id=temp[0] #AVANT : id = première colonne - id se réinitialise à quaque tour pour chaque ligne
                        res=[""]
                        res[0]=temp
                        if id not in reads.keys() : #If read's ID not in dictionary's keys
                            reads[id]=res      
                            compt1=compt1+1
                        else :
                            reads[id].append(temp)       #Else ID already on, .append !!!!!! on ajoute
                            compt1=compt1+1
                    else:
                        print("This is the header.") 
                file.close()     
          
#-----------------------------USE FOUNCTIONS TO LIMITE INDENTATION

#Creation of dictionary with flag's signification                  
defFlag={2048:"Supplementary alignment", 1024:"Read is PCR or optical duplicate", 512:"Read fails platform/vendor quality checks", 256:"Not primary alignment", 128:"Second in pair",64:"First in pair",32:"Mate reverse strand",16:"Read reverse strand", 8:"Mate unmapped", 4:"Read unmapped", 2:"Read mapped in proper pair", 1:"Read paired"}

#Here : method to know if a read is mapping or not. Return true or false. In same time count reads mapped.
def lectureFlag (dico) :
    compt2=0   #Counter for reads mapped
    compt3=0   #Counter for reads proped mapped
    idMapped=[]   #List mapped reads
    idProperMapped=[]    #List proper mapped list
    paire=0 #permet savoir si id a un ou deux reads
    for key in dico :   #For all our results      
        descr1=[]
        descr2=[]
        val1=int(dico[key][0][1])        #Flag value on SAM file
        idtemp1=dico[key][0][0]+"/1"
        if len(dico[key])>1 :
            val2=int(dico[key][1][1])
            idtemp2=dico[key][1][0]+"/2"
            paire=1
        else :
            paire=0
            val2=-1
        for key in defFlag :  #For each flag description
            if key<=val1 :       #If flag value it's less than flag or result of soustraction (to the secound run)
                descr1.append(defFlag[key])   #List of reads description 
                val1=val1-key      #Update the value to compare flag's description
            if paire==1 :
                if key<=val2 :
                    descr2.append(defFlag[key])
                    val2=val2-key
        if "Read unmapped" not in descr1 and len(descr1)!=0 : #Une fois que la foucle for juste au dessus est finie
            compt2=compt2+1
            idMapped.append(idtemp1)      #Stock ID of mapped reads     
        if "Read mapped in proper pair" in descr1 :
            compt3=compt3+1                       # Ici utiliser une troisième va pour voir si 
            idProperMapped.append(idtemp1.split("/")[0])      #Stock ID of proper mapped reads
        if paire==1 :
            if "Read unmapped" not in descr2 and len(descr2)!=0 :
                compt2=compt2+1
                idMapped.append(idtemp2)      #Stock ID of mapped reads         
        descr1=[] #Pour initialiser entre deux boucles
    return compt2, compt3, idMapped, idProperMapped
     




#Here : Change minimum limit to consider a read has a qood quality 
def choixSeuil () :
    essai1=0
    essai4=0
    val=30    #Initialisation of minimal quality limit by default
    while essai1!=1 :
        rep=input("La qualité de l'alignement d'un read se mesure entre 0 et 60. Le seuil minimum par défaut est de 30, où en général on considère la qualité comme bonne. /n Voulez vous la modifier ? Taper y ou n.") #WARNING :/n pour retour a la ligne. ok ???, + inclure un cf dicument readme pour expliquer calcul qualité ? idée : donner taux erreur relatif à 30 et calculer celui valeur donnée par user ?
        if rep!="y" and rep!="n" and rep!="yes" and rep!="no" :
            print("Votre réponse n'est pas celle attendue. Veillez taper y(yes) ou n (no), sans majuscule, svp.") 
            rep=input("Voulez-vous modifier la valeur par défaut du seuil minimum de qualité ?")
        else :
            essai1=1
            if rep=="y" or rep=="yes" :  
                while essai4!=1 :
                    val=int(input("Veuillez entrer la valeur de seuil désiré pour la qualité des reads sélectionnés."))
                    if val<0 or val>60 :
                        print("Le seuil doit être compris entre 0 et 60.")
                    else :
                        essai4=1
            else :
                print("OK. Les reads affichés auront une qualité d'au moins 30 sur 60.")
    return val


#Here : Count read mapped with the minimum quality and if there are tot mapped
def MappedQualOrTot(liste) :
    readMappeQual=[]
    readTotMappedQual=[]          
    prevI=''
    tailleRead=0
    for i in liste :
        idtemp1=i.split("/")[0]
        if idtemp1 in reads.keys() :
            if idtemp1!=prevI : #If it's first read of the pair
                if int(reads[idtemp1][0][4])>=seuil :               
                    readMappeQual.append(idtemp1+"/1") 
                    tailleRead=len(reads[idtemp1][0][9])
                    if reads[idtemp1][0][5]==str(tailleRead) +"M":
                        readTotMappedQual.append(reads[idtemp1][0][0]+"/1")
            if idtemp1==prevI : #If it's the second pair of the read
                if int(reads[idtemp1][1][4])>=seuil :
                    readMappeQual.append(idtemp1+"/2")
                    tailleRead=len(reads[idtemp1][1][9])
                    if reads[idtemp1][1][5]==str(tailleRead) +"M":
                        readTotMappedQual.append(reads[idtemp1][1][0]+"/2")
        prevI=idtemp1
    return readMappeQual, readTotMappedQual

#Here : Count read totaly mapped  (not with minimum qual required)
def compteur5 (dico):
    tailleRead=0
    idTotMapped=[]
    for key in dico :
        tailleRead=len(dico[key][0][9])
        if dico[key][0][5]==str(tailleRead) +"M":
            idTotMapped.append(dico[key][0][0]+"/1")
    if len(dico[key])>1 :   #If there is a paired read
        tailleRead=len(dico[key][1][9])
        if dico[key][1][5]==str(tailleRead) +"M":
            idTotMapped.append(dico[key][1][0]+"/2")
    return idTotMapped


#Here : Ask user if he want export result in new text file
def exportReads():
    essai2=0
    essai3=0
    while essai2!=1 :
        aff=input("Do you want export this file analysing results on a new text file ?")
        if aff!="y" and aff!="n" and aff!="yes" and aff!="no" :
            print("Votre réponse n'est pas celle attendue. Veillez taper y(yes) ou n (no), sans majuscule, svp.") 
            aff=input("Tapez votre réponse en minuscule.")
        else :
            essai2=1
            if aff=="y" or aff=="yes" :
                nameFile=input("Write the text file's name wanted please.")
                while essai3!=1 :
                    if "." in nameFile :
                        print("File's name mustn't containt the . symbol (don't try to put file's extension !).")
                        nameFile=input("Write the text file's name wanted please.")
                    else :
                        essai3=1
                f=open(nameFile+".txt","w")
                f.write("Résultats de l'analyse du fichier "+sys.argv[1].split("/")[-1]+" par le programme recapSAM.py")
                f.write("\n \n recap.SAM count "+str(compt1)+" reads on the SAM file.")
                f.write("\n \n On this reads there are "+str(len(readsMapped))+" reads are mapped. These represente "+str(quantMappedOnTot)+"% of all reads. On this mapped reads there are :")
                f.write("\n - "+str(len(readsProperMapped))+" pairs of proper mapped reads")
                f.write("\n - "+str(len(readsMappedQual))+" mapped with a minimal quality of "+str(seuil)+". That represente "+str(quantMappedQualOnTot)+"% of total reads and "+str(quantMappedQualOnMapped)+"% of mapped reads.")
                f.write("\n - "+str(len(readsTotMapped))+" reads totally mapped (not difference between read's sequence and the reference sequence on alignment).")
                if len(readsTotMappedQual)!=0 :
                    f.write("\n \n On the mapped reads with the minimum quality required there are "+str(len(readsTotMappedQual))+" reads totally mapped.")
                    f.write("\n That represente : \n- "+str(quantTotMappQualOnReads)+"% of all reads, \n- "+str(quantTotMappQualOnMapp)+"% of mapped reads with the minimal required quality, \n- "+str(quantTotMappQualOnTotMapp)+"% of totally mapped reads (there are "+str(len(readsTotMapped))+" reads totally mapped on this file).")
                else :
                    f.write("\n \n There are not reads totally mapped with the minimum quality "+str(seuil)+" in this file.")
                f.close()  
                
           


#-----------------------------CORPS DU PROGRAMME

#test ici pour moins d'indentation 
if any(reads) : #If file can was readed with success : because any vérifie que reads ne soit pas vide
    print("Le nombre de reads totals dans votre fichier est de "+ str(compt1))
    count2, count3, readsMapped, readsProperMapped=lectureFlag(reads)  #à gauche les variables ou sont stocké ce que la fonction retourne
    quantMappedOnTot="{:.2f}".format((count2/compt1)*100)
    print("Number of mapped reads on these reads : "+str(count2)+".  That represente "+str(quantMappedOnTot)+"% of reads.") #nbMappe is result of  our second counter
    print("Number of PROPER mapped pairs of reads on this mapped reads is : "+str(count3))
    seuil=choixSeuil()   #User choice value of minimum quality he wants
    print("This value of minimal quality mean the highest probability of the mapping position is wrong is equal to "+str(10**(-seuil/10)))
    readsMappedQual, readsTotMappedQual=MappedQualOrTot(readsMapped)
    quantMappedQualOnMapped="{:.2f}".format((len(readsMappedQual)/count2)*100)
    quantMappedQualOnTot="{:.2f}".format((len(readsMappedQual)/compt1)*100)
    print("Les reads mappés d'une qualité minimale de "+str(seuil)+" sont au nombre de :"+str(len(readsMappedQual))) 
    print("That represente "+str(quantMappedQualOnTot)+"% of total reads and "+str(quantMappedQualOnMapped)+"% of mapped reads.")
    readsTotMapped=compteur5(reads)   
    if len(readsTotMappedQual)!=0 :
        quantTotMappQualOnReads="{:.2f}".format((len(readsTotMappedQual)/compt1)*100)
        quantTotMappQualOnMapp="{:.2f}".format((len(readsTotMappedQual)/len(readsMappedQual))*100)
        quantTotMappQualOnTotMapp="{:.2f}".format((len(readsTotMappedQual)/len(readsTotMapped))*100)
        print("Number of reads totaly mapped with the minimum quality required ("+str(seuil)+") is : "+str(len(readsTotMappedQual)))
        print("These reads mapped with minimal quality represented "+str(quantTotMappQualOnReads)+" % of all reads on file, "+str(quantTotMappQualOnReads)+" % of reads mapped with minimum quality required, and "+str(quantTotMappQualOnTotMapp)+" % of reads totaly mapped.")
    else :
        print("There are not reads totally mapped with the minimum quality "+str(seuil)+" in this file.")
    exportReads()
else :
    print("Votre entrée ne peut être traitée. Veuiller entre en premier argument (après le nom du programme) un fichier SAM non-vide.")


















