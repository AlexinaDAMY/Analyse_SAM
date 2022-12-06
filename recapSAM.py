#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Import necessary module
import sys, os

#-----------------------------ANALIZED FILE : SAM file ==> READING

reads={}   
compt1=0   
if len(sys.argv)>1 :
    if os.path.isdir(sys.argv[1]) :
        print("This program can't process many files on a directory. Please input only one file.")
    else :                                                                
        entree=sys.argv[1]
        fichierEntre=entree.split("/")[-1]       
        extFichier=fichierEntre.split(".")[-1]   
        if extFichier!="sam" and extFichier!="SAM" : 
            print("This program it's designed for SAM file.")
        else :                                                           
            if os.path.getsize(sys.argv[1])==0 :     
                print("Analysing file is empty. The program recapSAM.py can't run.")
            else :
                file=open(sys.argv[1],"r")
                lines=file.readlines() 
                for line in lines :
                    temp=line.split("\t") 
                    if temp[0][0]!='@' : 
                        id=temp[0] 
                        res=[""]
                        res[0]=temp
                        if id not in reads.keys() : 
                            reads[id]=res      
                            compt1=compt1+1
                        else :
                            reads[id].append(temp)      
                            compt1=compt1+1
                    else:
                        print("This is the header.") 
                file.close()     
          
#-----------------------------USE FOUNCTIONS TO LIMITE INDENTATION

#Creation of dictionary with flag's signification                  
defFlag={2048:"Supplementary alignment", 1024:"Read is PCR or optical duplicate", 512:"Read fails platform/vendor quality checks", 256:"Not primary alignment", 128:"Second in pair",64:"First in pair",32:"Mate reverse strand",16:"Read reverse strand", 8:"Mate unmapped", 4:"Read unmapped", 2:"Read mapped in proper pair", 1:"Read paired"}

#Here : method to know if a read is mapping or not and if it's mapped in proper pair. 
def lectureFlag (dico) :
    compt2=0  
    compt3=0   
    idMapped=[]   
    idProperMapped=[]    
    paire=0 
    for key in dico :       
        descr1=[]
        descr2=[]
        val1=int(dico[key][0][1])        
        idtemp1=dico[key][0][0]+"/1"
        if len(dico[key])>1 :
            val2=int(dico[key][1][1])
            idtemp2=dico[key][1][0]+"/2"
            paire=1
        else :
            paire=0
            val2=-1
        for key in defFlag :  
            if key<=val1 :      
                descr1.append(defFlag[key])   
                val1=val1-key      
            if paire==1 :
                if key<=val2 :
                    descr2.append(defFlag[key])
                    val2=val2-key
        if "Read unmapped" not in descr1 and len(descr1)!=0 : 
            compt2=compt2+1
            idMapped.append(idtemp1)           
        if "Read mapped in proper pair" in descr1 :
            compt3=compt3+1                   
            idProperMapped.append(idtemp1.split("/")[0])     
        if paire==1 :
            if "Read unmapped" not in descr2 and len(descr2)!=0 :
                compt2=compt2+1
                idMapped.append(idtemp2)            
        descr1=[] 
    return compt2, compt3, idMapped, idProperMapped
     




#Here : Change minimum limit to consider a read has a qood quality 
def choixSeuil () :
    essai1=0
    essai4=0
    val=30    
    while essai1!=1 :
        rep=input("Alignement quality of one read is included between 0 and 60. There the minimum quality value by default is 30 (in general it's considerated like a good value). /n Do you want change this value ? Write y or n.")
        if rep!="y" and rep!="n" and rep!="yes" and rep!="no" :
            print("Your answer it's not that expected. Please write y (yes) or n (no), without capital letter.") 
            rep=input("Do you want change this value ?")
        else :
            essai1=1
            if rep=="y" or rep=="yes" :  
                while essai4!=1 :
                    val=int(input("Please write your value for minimal quality of mapped reads selected."))
                    if val<0 or val>60 :
                        print("This value must be included between 0 and 60 !")
                    else :
                        essai4=1
            else :
                print("Okay. I select reads with a minimal quality equal to 30.")
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
            if idtemp1!=prevI : 
                if int(reads[idtemp1][0][4])>=seuil :               
                    readMappeQual.append(idtemp1+"/1") 
                    tailleRead=len(reads[idtemp1][0][9])
                    if reads[idtemp1][0][5]==str(tailleRead) +"M":
                        readTotMappedQual.append(reads[idtemp1][0][0]+"/1")
            if idtemp1==prevI : 
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
        if len(dico[key])>1 :   
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
            print("Your answer it's not that expected. Please write y (yes) or n (no), without capital letter.") 
            aff=input("Do you want export results on text file ?")
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
                f.write("Analysing SAM file results "+sys.argv[1].split("/")[-1]+" by the program recapSAM.py")
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
                
           


#-----------------------------PROGRAM'S BODY
 
if any(reads) : 
    print("Number of reads on your file it's "+ str(compt1))
    count2, count3, readsMapped, readsProperMapped=lectureFlag(reads)  
    quantMappedOnTot="{:.2f}".format((count2/compt1)*100)
    print("Number of mapped reads on these reads : "+str(count2)+".  That represente "+str(quantMappedOnTot)+"% of reads.") 
    print("Number of PROPER mapped pairs of reads on this mapped reads is : "+str(count3))
    seuil=choixSeuil() 
    print("This value of minimal quality mean the highest probability of the mapping position is wrong is equal to "+str(10**(-seuil/10)))
    readsMappedQual, readsTotMappedQual=MappedQualOrTot(readsMapped)
    quantMappedQualOnMapped="{:.2f}".format((len(readsMappedQual)/count2)*100)
    quantMappedQualOnTot="{:.2f}".format((len(readsMappedQual)/compt1)*100)
    print("Number of mapped reads with one minimal quality of "+str(seuil)+" is :"+str(len(readsMappedQual))) 
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
        print("But there are "+str(len(readsTotMapped))+" reads totaly mapped with a lower quality.")
    exportReads()
else :
    print("Your entry can't be process. Please write on first argument (after the program name) one SAM file path, with file not empty.")


















