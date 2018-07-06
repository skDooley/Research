'''
 Tools to find, annotate, and annalyze CRISPR Loci
'''
#import sys
import tempfile
from Bio import pairwise2
from Bio.SeqIO import write as swrite
from Bio.SubsMat import MatrixInfo as matlist
from easyFunctions import *
#from Bio.Blast.Applications import NcbiblastnCommandline
from GetOrfs import *
#from Bio.SeqIO import write as writeFasta
#from Bio.SeqIO import write
#from Bio.Alphabet import generic_dna

class FileWrapper:
    def __init__(self,filename):
        self.index = -1
        self.lines = open(filename,"r").readlines()
    def readline(self):
        self.index += 1
        return self.lines[self.index]
    def has_more_lines(self): return len(self.lines)-1 != (self.index) 
    
matrix = matlist.blosum62
gap_open = -10
gap_extend = -.8

def Orientation(coordinate): return coordinate.start > coordinate.end
def CalcPercIdent(seqA,seqB):
    seq1,seq2 = alignSequence(seqA,seqB)
    matches = 0
    for index in range(len(seq1)):
        if seq1[index] == seq2[index]:matches += 1
    return float(matches)/len(seqA)*100

def alignSequence(refSeq,fragment):
    probe = refSeq #scrubSequence(refSeq)
    ext_plus_lig = fragment
    alns = pairwise2.align.globalds(probe, ext_plus_lig, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_probe, aln_arms, score, begin, end = top_aln
    return [aln_probe, aln_arms]

class repeatStruct:
    def __init__(self,id,seq): self.id = id; self.seq = seq
    def __hash__(self): return hash(str(self.id)+self.seq)
    def __eq__(self,other): return self.id == other.id and self.seq == other.seq

def findNCoords(seq):
    seq = seq.upper()
    index = -1
    coords = []
    while seq.find("N"*20,index+1) > -1:
        index = seq.find("N"*20,index+1)
        start_index = index
        while index<len(seq) and seq[index] == 'N': index+=1
        coords.append(Coordinate(start_index,index))
#    if len(coords)>0:
#        print "\tThese are the ncoords",
#        for cord in coords: print "\t\t",cord
    return coords

class CRISPRLocus:
    def __init__(self,type):
        self.name = ""
        self.assembly_size = 0
        self.CRISPRarrayLocations = []
        self.repeats = {}
        self.consensusRepeat = []
        self.repeatCoords = []
        self.spacers = set()
        self.spacerCoords = []
        self.antiRepeats = set()
        self.possibleTracrs = set()
        self.type = type
        self.orfs = None
        self.min,self.max = -1,-1
        self.hasNRegion = False
    
    def annotate(self,seq):
        minCoord, maxCoord = len(seq), 0
        
        #BLAST for possible protein matches in NR
        seqFile = open("tmp/genomicSequence.fasta","w")
        seqFile.write(">%s\n%s\n" % (self.name,seq))
        seqFile.close()
        orfs = getOrfs("tmp/genomicSequence.fasta",90)
        blastResults = "data/Annotations/blast/%s_orfs.blastout" % (self.name)
        fileHandle = FileWrapper(blastResults)
        fileHandle.readline()
        while fileHandle.has_more_lines():
            line = fileHandle.readline().strip()    
            if "# Query: " in line: 
                ID = line[len("# Query: "):line.find("|")] #Query
                fileHandle.readline() #database
                line = fileHandle.readline().strip() #hits found or fields if found something
                if line == "# 0 hits found":continue
                fileHandle.readline()
                rec = fileHandle.readline().strip().split("\t")
                ORF_ID,sciName,description = rec[0],rec[-2],rec[-1]
                if "hypothetical" in description.lower():continue
                words = description.lower().split(" ")
                goodWords = set(["crispr","endonuclease","cas9","csn1","cas1","cas2","crrna","tracrrna","crispr-associated","csn2","dna2/cas4","cas4","cas1-like", "cas1/cas4","precrrna", "crrna", "tracr", "tracrrna", "tnpb", "pre-crrna"])
                goodWordFound = False
                for word in words:
                    goodWordFound = (word in goodWords)
                    if goodWordFound:break
                if not goodWordFound:continue
                found = False
                for feature in orfs.records[ID].features:
                    if len(feature.qualifiers['label']) ==0:continue
                    if feature.qualifiers['label'][0] == ORF_ID and not found:
                        feature.qualifiers['label'] = description[description.find(" ")+1:description.find("[")]
                        feature.qualifiers['note'] = description+" "+ORF_ID
                        feature.qualifiers['blast_hit'] = "True"
                        minCoord =  min(feature.location.start,feature.location.end,minCoord)
                        maxCoord =  max(feature.location.start,feature.location.end,maxCoord)
                        found = True
                        break
        for id in orfs.records:
            remove = []
            for i in range(len(orfs.records[id].features)):
                feature = orfs.records[id].features[i]
                if 'blast_hit' not in feature.qualifiers:remove.append(i)
            for index in remove[::-1]: del orfs.records[id].features[index]
        
        #ADD the CRISPR ARRAY HITs/BLAST
        repeatIndex = 0
        for repeat in self.repeatCoords:
            fid = "Repeat_%i" % (repeatIndex)
            minCoord =  min(repeat.start,repeat.end,minCoord)
            maxCoord =  max(repeat.start,repeat.end,maxCoord)
            repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="repeat", strand=0)
            repeatFeature.qualifiers['label'] = [fid]
            orfs.records[self.name].features.append(repeatFeature)
            repeatIndex+=1
        repeatIndex = 0
        for repeat in self.possibleTracrs:
            fid = "Theoretical antirepeat_%i" % (repeatIndex)
            minCoord =  min(repeat.start,repeat.end,minCoord)
            maxCoord =  max(repeat.start,repeat.end,maxCoord)
            repeatFeature = SeqFeature(FeatureLocation(min(repeat.start,repeat.end), max(repeat.start,repeat.end)), type="terminator", strand=0)
            repeatFeature.qualifiers['label'] = [fid]
            orfs.records[self.name].features.append(repeatFeature)
            repeatIndex+=1
        
        #Get N coordinates
        NCoords = findNCoords(seq)
        for coord in NCoords:
            if ((minCoord-500) <= coord.start and coord.start <= (maxCoord+500)) or (max(minCoord-500,0) <= coord.end and coord.end <= (maxCoord+500)):
                minCoord =  min(coord.start,coord.end,minCoord)
                maxCoord =  max(coord.start,coord.end,maxCoord)
                ambigSeq = SeqFeature(FeatureLocation(coord.start,coord.end), type="ambigious_seq", strand=0)
                ambigSeq.qualifiers['label'] = ["AmbigiousSequence"]
                orfs.records[self.name].features.append(ambigSeq)
                self.hasNRegion = True
        
        self.min,self.max, self.orfs = minCoord, maxCoord,orfs
        

    def confirmedAnti(self,anti_repeat):
        remove = set()
        foundAnti = False
        for rep in self.antiRepeats:
            foundAnti = anti_repeat.overlaps(rep)
            if not foundAnti:remove.add(rep)
            if foundAnti: break
        if foundAnti: return foundAnti
#        for coord in remove: self.antiRepeats.remove(coord)
        return False
    
    def combineResults(self,other):
        if self.type != "Piler": self.consensusRepeat = other.consensusRepeat[0]
        self.type = "Combined"
        #Compare CRISPR Array Locations
        locA_index = 0
        for locationA in other.CRISPRarrayLocations:
            foundLocation = False
            locB_index = 0
            for locationB in self.CRISPRarrayLocations:
                foundLocation = (locationA == locationB)
                if foundLocation:
                    self.CRISPRarrayLocations[locA_index].start = min(self.CRISPRarrayLocations[locA_index].start,self.CRISPRarrayLocations[locB_index].start)
                    self.CRISPRarrayLocations[locA_index].end = max(self.CRISPRarrayLocations[locA_index].end,self.CRISPRarrayLocations[locB_index].end)
                    break
                locB_index += 1
            locA_index += 1
            if not foundLocation: self.CRISPRarrayLocations.append(locationA)
        
        #Combine anti-repeats
        for anti_repeatA in other.antiRepeats:
            #If the anti-repeat is in the crispr array of the other found array then continue on to the next repeat
            inCRISPRArray = False
            for location in self.CRISPRarrayLocations:
                inCRISPRArray = (location.contains(anti_repeatA) or anti_repeatA.overlaps(location))
                if inCRISPRArray:break
            if inCRISPRArray:continue
            
            #Check to see if the anti-repeat is a duplicate of an existing anti    
            foundAntiRepeat = False
            for anti_repeatB in self.antiRepeats:
                foundAntiRepeat = (anti_repeatA == anti_repeatB)
                if foundAntiRepeat: break
            if not foundAntiRepeat: self.antiRepeats.add(anti_repeatA)
        
        #Combine Repeats
        lenRepCoords = len(self.repeatCoords)
        if lenRepCoords != len(set(self.repeatCoords)):
            print lenRepCoords,len(set(self.repeatCoords))
            raise Exception("How can a repeat coord have a start position that is the same?")
        theRepeats = set(self.repeatCoords)
        numFound = 0
        for coordA in other.repeatCoords:
            if coordA in theRepeats:pass
            else:
                foundCoord = False
                for coordB in theRepeats:
                    foundCoord = (coordA == coordB)
                    if foundCoord:break
                if not foundCoord: self.repeats["NewRepeat_%i" %(numFound)] = coordA
    
    def checkExtendsArray(self,coord,crRNALen):
        print "\t\t",coord,
        if self.CRISPRarrayLocations[0].distance(coord) <= max(crRNALen,100):
            self.extendArray(coord)
            print "Extends" #"\tExtends = TRUE for ",coord
        else:
            print "False"
            return coord
    
    def checkArrayBoundaries(self, blast_results):
        print self.name + " CRISPR Array:",self.CRISPRarrayLocations[0]
        before,after = [],[]
        for hit in blast_results:
            if hit.start < self.CRISPRarrayLocations[0].start: before.append(hit)
            else:after.append(hit)  

        retHits = []
        avgSpacerLen = 0
        for spacer in self.spacers: avgSpacerLen += len(spacer)
        avgSpacerLen = int(avgSpacerLen/float(len(self.spacers))) 
        crRNALen = avgSpacerLen + len(self.consensusRepeat[0])
        print "\tAverage spacer len: %i + Repeat len: %i = %i basepairs" % (avgSpacerLen,len(self.consensusRepeat[0]),crRNALen)
        #For the each of the blast hits
        ##Ignore any hits that are completely contained within the crispr array
        ##If there is overlap but not complete containment, update the array
        ##If the coordinate doesn't overlap, check to see if it extends
        ###If it extends update the array and the repeats
        ###Else add to list of hits to return
        print "\tBefore:"
        for coord in before[::-1]: #For the each of the spacially ORDERED blast hits
            if self.CRISPRarrayLocations[0].contains(coord):print "\t",coord,"Within"; continue ##Ignore any hits that are completely contained within the crispr array
            if self.CRISPRarrayLocations[0].overlaps(coord): 
                print "\t\t",coord,"Extends"
                coord = self.extendArray(coord) #Check for a blast hit outside of the crispr array 
            else: coord = self.checkExtendsArray(coord,crRNALen) #Check that the coordinate is completly contained within the CRISPR
            if coord is not None:
                #print "Before", coord,self.CRISPRarrayLocations[0].contains(coord),self.CRISPRarrayLocations[0].overlaps(coord)
                retHits.append(coord)

        print "\tAfter:"
        for coord in after: 
            if self.CRISPRarrayLocations[0].contains(coord):print "\t\t",coord,"Within"; continue ##Ignore any hits that are completely contained within the crispr array
            if self.CRISPRarrayLocations[0].overlaps(coord): 
                print "\t\t",coord,"Extends"
                coord = self.extendArray(coord) #Check for a blast hit outside of the crispr array 
            else: coord = self.checkExtendsArray(coord,crRNALen) #Check that the coordinate is completly contained within the CRISPR
            if coord is not None:
                retHits.append(coord)
        print
        print "\tUpdated CRISPR Array:",self.CRISPRarrayLocations[0]
        print "\tPossible TRACRS:"
        print "\t\t",self.repeatCoords[0]
        print "\t\t",self.repeatCoords[-1]
        for hit in retHits: print "\t\t",hit
        print
        self.possibleTracrs.add(self.repeatCoords[0])
        self.possibleTracrs.add(self.repeatCoords[-1])
        for coord in retHits:self.possibleTracrs.add(coord)
        return retHits, crRNALen    
        
    def getTerminusRepeats(self,seq,index):
        index+=1
        retStr = ""
        fastaSeq = ">Seq_%s_%s_%i_%i_%i\n%s\n" 
        startCoord = self.repeatCoords[0]
        endCoord = self.repeatCoords[-1]
        retStr += fastaSeq % ("E",self.consensusRepeat[0],min(startCoord.start,startCoord.end),max(startCoord.end,startCoord.start),index,seq[min(startCoord.start,startCoord.end)-500:max(startCoord.end,startCoord.start)])
        index+=1
        retStr += fastaSeq % ("S",self.consensusRepeat[0],min(startCoord.start,startCoord.end),max(startCoord.end,startCoord.start),index,seq[min(endCoord.start,endCoord.end):max(endCoord.end,endCoord.start)+500])
        return retStr,index
            
    def extendArray(self,newCoord):
        if newCoord.start < self.CRISPRarrayLocations[0].start:
            self.CRISPRarrayLocations[0].start = newCoord.start
            self.repeatCoords = [newCoord] + self.repeatCoords
        else: 
            self.CRISPRarrayLocations[0].end = newCoord.end
            self.repeatCoords.append(newCoord)

    def addrepeat(self,repeat):
        try: self.repeats[len(self.CRISPRarrayLocations)][repeat]+=1
        except: 
            try: self.repeats[len(self.CRISPRarrayLocations)][repeat] = 1
            except: self.repeats[len(self.CRISPRarrayLocations)] = {repeat:1}
        
    def add_crRNA(self,pos,repeat,spacer):
        self.spacers.add(spacer)
        self.addrepeat(repeat)
        self.repeatCoords.append(Coordinate(int(pos)-1,int(pos) + len(repeat)-1))
        self.spacerCoords.append(Coordinate(int(pos) + len(repeat),int(pos) + len(repeat) + len(spacer)-1))
    
    def combineRepeats(self,repeatA,repeatB):
        repeatAID,repeatASeq,repeatBID,repeatBSeq,countA,countB = "","","","",0,0
        for ID,seqDict in self.repeats.iteritems():
            for repeat in seqDict:
                if repeat == repeatA: 
                    repeatAID = ID
                    repeatASeq = repeat
                    countA = self.repeats[ID][repeat]
                    break
                elif repeat == repeatB: 
                    repeatBID = ID
                    repeatBSeq = repeat
                    countB = self.repeats[ID][repeat]
                    break
        if repeatASeq != "" and repeatBSeq != "":
            if countA > countB:
                self.repeats[repeatAID][repeatASeq] += countB
                del self.repeats[repeatBID][repeatBSeq]
                if len(self.repeats[repeatBID]) == 0: del self.repeats[repeatBID]
            else:
                self.repeats[repeatBID][repeatBSeq] += countA
                del self.repeats[repeatAID][repeatASeq]
                if len(self.repeats[repeatAID]) == 0: del self.repeats[repeatAID]
    
    def removeDuplicates(self):
        remove = set()
        for IDA,seqDictA in self.repeats.iteritems():
            for IDB,seqDictB in self.repeats.iteritems():
                if IDA != IDB:
                    for seqA in seqDictA:
                        for seqB in seqDictB:
                            if seqA == seqB: 
                                if self.repeats[IDA][seqA] > self.repeats[IDB][seqB]: remove.add(repeatStruct(IDB,seqB))
                                else: remove.add(repeatStruct(IDA,seqA))
        if len(remove) != 0:
            for coord in remove:
                del self.repeats[coord.id][coord.seq]
                if len(self.repeats[coord.id]) == 0: del self.repeats[coord.id]
        
    def calc_consensus(self): 
        maxCount, maxSeq = 0, ""
        for id, seqDict in self.repeats.iteritems():
            for sequence, count in seqDict.iteritems(): 
                if count>maxCount: 
                    maxSeq = sequence
                    maxCount = count
#                print sequence, count
#                allRepeats.add(sequence)
        self.consensusRepeat=[maxSeq]
#        self.consensusRepeat=list(allRepeats)
#        #Strategy 1 (Depricated): 
#        '''
#        If there is is only one crispr array location, then compute a consensus
#        by aligning all of the repeats and taking the IUPAC for the allelels, ignore
#        gaps if  gap_percentage >= .75 or gap_percentage <= .25 (typical cases meaning really rare or almost completly gapped
#        '''
#        if self.type == "Piler": raise Exception("This should not happen")
#        if len(self.CRISPRarrayLocations) == 1:
#            alignments = alignSeqs(self.repeats)
#            seqs=[]
#            for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#            self.consensusRepeat.append(calc_consensus_seq(seqs,True))
#        else:
#            for id,seqs in self.repeats.iteritems():
#                alignments = alignSeqs(seqs)
#                seqs=[]
#                for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#                self.consensusRepeat.append(calc_consensus_seq(seqs,True))
#                
#            AllPossibleRepeats = set()
#            for id, repeatDict in self.repeats.iteritems():
#                for repeat in repeatDict: AllPossibleRepeats.add(repeat)
#            ignore = set()
#            for repeatA in AllPossibleRepeats:
#                for repeatB in AllPossibleRepeats:
#                    if repeatA != repeatB and ("%s_%s" % (repeatA,repeatB) not in ignore):
#                        percIdent = CalcPercIdent(repeatA,repeatB)
#                        if percIdent > 90.0: self.combineRepeats(repeatA,repeatB)
#                        ignore.add("%s_%s" % (repeatA,repeatB))
#                        ignore.add("%s_%s" % (repeatB,repeatA))
#            self.removeDuplicates()
#            if len(self.CRISPRarrayLocations) == 1 and len(self.repeats[1]) == 1: #expect to blow up when dupes in 1-2 and 1 is deleted cause more in 2
#                for id,seqDict in self.repeats.iteritems():
#                    self.consensusRepeat = seqDict.keys()
#                    if len(self.consensusRepeat) != 1: raise Exception("Problem here " +str(self.repeats))
#                    return
#            alignments = alignSeqs(self.repeats)
#            for alnmnt in alignments: seqs.append(str(alnmnt.seq))
#            for id,seqDict in self.repeats.iteritems():
#                for seq in seqDict:
#                    self.consensusRepeat.append(seq)
#            self.consensusRepeat = calc_consensus_seq(seqs)
            

    def resolveRepeats(self,summary_line):
        '''
        This method is to resolve repeats from pilercr which are given as mathes
        being '.' and SNPs given as the actual base. @parameter summary_line looks like:
                29      46              29                GTTGTGAATAGCTTTCAAAATTGTATCTTAGTAGATGATTCACAGG
        NumRepeats RepeatLen   %id(blank)  SpacerLen    Consensus
        '''
        self.consensusRepeat.append(summary_line.split("    ")[-1])
        #TODO: Could resolve repeats into sequence but for now will leave them as ....A....T......

    def findAntiRepeats(self,record):
        if len(self.consensusRepeat) == 0 and len(self.repeatCoords)!= 0: self.calc_consensus()
        if len(self.consensusRepeat) == 0: raise Exception("No repeats for", self.name, len(self.repeatCoords))
        swrite([record],"tmp/assemblyFasta.fa","fasta")
        repeatFasta = MakeFasta(tempfile.NamedTemporaryFile(delete=False))
        counter = 0
        for repeat in self.consensusRepeat:
            repeatFasta.write("%s_%i" %(self.name,counter), repeat)
            counter += 1
        repeatFasta.close()
        results = BLAST(repeatFasta.name,"tmp/assemblyFasta.fa","tmp/blast_results.out")
        repeatFasta.die()
        remainingResults = self.compareResults(results) # This assumes that a hit is the same if it is contained or within the CRISPR
        for result in remainingResults: self.antiRepeats.add(result)
            #if Orientation(result): self.antiRepeats.append(result)
            #else: self.repeatCoords.append(result)

    def compareResults(self,results):
        #Assumption: The anti repeat is not in the crispr array and doesn't overlap it
        retList = []
        for result in results:
            foundInArray = False
            overlapsLocation = False
            for location in self.CRISPRarrayLocations:
                foundInArray = location.contains(result)
                if foundInArray: break
                overlapsLocation = location.overlaps(result)
                if overlapsLocation: break
            minSeparationFromArray = min(abs(location.start-result.start),abs(location.start-result.end),abs(location.end-result.end),abs(location.end-result.start))
            if not foundInArray and minSeparationFromArray > 1 and not overlapsLocation: retList.append(result)
        return retList

    def consensus(self):return self.calc_consensus()
    def __str__(self): 
#        crRNAs = ""
#        for index in range(0,len(self.repeatCoords)): crRNAs += "%i\t%i\t%i\t%i\n" %(self.repeatCoords[index].start, self.repeatCoords[index].end, self.spacerCoords[index].start,self.spacerCoords[index].end)
        #AssemblyName \t CRISPR_Start \t CRISPR_End \t Consensus Repeat \n
#        print "%s\t%i\t%i\t%s\n" % (self.name,self.CRISPRarrayLocations.start,self.CRISPRarrayLocations.end,self.consensusRepeat)
        retstr = ""
#        if type(self.consensusRepeat) == type([]):self.consensusRepeat = self.consensusRepeat[0]
        for crArray in self.CRISPRarrayLocations: retstr += "%s\t%i\t%i\t%s\n" % (self.name,crArray.start,crArray.end,str(self.consensusRepeat))         
        return retstr

from pandas import DataFrame
class TRACR_Locus(dict):
    def __init__(self, name): 
        self.name = name
        self.structuralHit = {} # TRACRStructureName:[Coordinates]
        self.overlaps = False
        self.topHit = None
    def _setTopHit(self):
        #Only 1 hit
        if len(self.structuralHit)==1:
            [(self.topHit, v)] = self.structuralHit.items()
            return
        #Check overlap
        df = DataFrame(index=self.structuralHit,columns=self.structuralHit)
        for hitName1,hit1 in self.structuralHit.iteritems():
            for hitName2,hit2 in self.structuralHit.iteritems():
                if hitName1 == hitName2: continue
                if hit1.location.overlaps(hit2.location):
                    dif = hit1.location - hit2.location
                    if df[hitName2][hitName1] != dif: df[hitName1][hitName2] = dif
        maxID1,maxID2 = df.max(axis=1).idxmax(),df.max(axis=0).idxmax()
        #If no overlap
        if maxID1 != maxID1:
            print df
            return #This case is interesting because overlap of blast hit with tracr but no overlap of tracrs
        #else if overlap
        if len(self.structuralHit[maxID1].location)>len(self.structuralHit[maxID2].location): self.topHit = maxID1
        else: self.topHit = maxID2    

    def structureHit(self,structureHit):
        for structList in structureHit:
            for structHit in structList:self.structuralHit[structHit.structName] = structHit
    def addBLASTHit(self,blastHit):
        for hitID,coord in blastHit: #hit = [StructName,coordinate]
            if hitID in self.structuralHit:
                if self.structuralHit[hitID].location.overlaps(coord): 
                    self.structuralHit[hitID].overlapsBLAST = True
                    self.overlaps = True
        self._setTopHit()
    def addHit(self,hitType,structureName,coordinate):
        try: self[hitType][structureName] = coordinate
        except: self[hitType] = {structureName:coordinate}
    def __hash__(self):return hash(self.name)

class FileWrapper:
    def __init__(self,filename):
        self.filename = filename
        self.index = -1
        self.lines = open(filename,"r").readlines()
    def readline(self):
        self.index += 1
        try:
            return self.lines[self.index]
        except:
            print "Failed at index", self.index, "for file", self.filename
            die
    def has_more_lines(self): return len(self.lines)-1 != (self.index)

class PilerCRReader(dict):
    def __init__(self, filename): 
        try: self.parse_file(FileWrapper(filename))
        except: pass
    
    def parse_file(self,file):
        for i in range(11):file.readline()
        while file.has_more_lines():
            newLocus = CRISPRLocus("Piler")
            line = file.readline()#Array # in file
            if "SUMMARY BY SIMILARITY" in line: break
            newLocus.name = file.readline().strip().replace(">","") #Contig ID
            if " " in newLocus.name: newLocus.name=newLocus.name[:newLocus.name.find(" ")]
            file.readline() #Blank Line
            file.readline() #Header
            file.readline() #Divider line made of =======
            line = file.readline() #1rst repeat
            while "==" not in line:
                crRNAInfo = line.strip().split("    ")
                if len(newLocus.CRISPRarrayLocations) == 0: newLocus.CRISPRarrayLocations.append(Coordinate(int(crRNAInfo[0]), int(crRNAInfo[0])))
                else: newLocus.CRISPRarrayLocations[0].end = int(crRNAInfo[0]) + len(crRNAInfo[3])
                newLocus.add_crRNA(int(crRNAInfo[0]), crRNAInfo[3], crRNAInfo[4])
                self[newLocus.name] = newLocus
                line = file.readline()
            line = file.readline().strip()#Summary Line
            newLocus.resolveRepeats(line)
            file.readline() #Blank Line
            file.readline() #Blank Line

class MinCEDReader(dict):
    def __init__(self, filename): self.parse_file(FileWrapper(filename))
    
    def readLocus(self, line, newLocus,file):
        while "--------" not in line:
            repeat = ""
            try: pos, empty, repeat, spacer, lengths = line.strip().split("\t")
            except ValueError: pos, empty, repeat = line.strip().split("\t")
            newLocus.add_crRNA(pos, repeat, spacer)
            line = file.readline()
                
    def parse_file(self,file):
        line = ""
        while file.has_more_lines():
            newLocus = CRISPRLocus("MinCED")
            if "Sequence " in line: newLocus.name, newLocus.assembly_size = self.getSeqName(line) #Line 1
            else: newLocus.name, newLocus.assembly_size = self.getSeqName(file.readline().strip()) #Line 1
            file.readline() #Line 2 is always blank
            coordLine = file.readline().strip().split(" ") #Line 3
            newLocus.CRISPRarrayLocations.append(Coordinate(int(coordLine[5]), int(coordLine[7])))
            file.readline() #Line 4 is a header line: POSITION REPEAT SPACER
            file.readline() #Line 5 is a divider line w/a bunch of ---
            line = file.readline() #ResultLines
            self.readLocus(line,newLocus,file)
            file.readline() #Summary line
            file.readline() #Blank line
            line = file.readline().strip() #line with time to find repeats
            while "Time to" not in line:
                coordLine = line.split(" ")
                newLocus.CRISPRarrayLocations.append(Coordinate(int(coordLine[5]), int(coordLine[7])))
                file.readline()
                file.readline()
                line = file.readline()
                self.readLocus(line,newLocus,file)
                file.readline()#Summary line
                file.readline()#Blank line
                line = file.readline()#line with time to find repeats
            #2 blank lines
            file.readline()
            file.readline()
            self[newLocus.name] = newLocus
    
    def calulate_consensus_seqs(self):
        for locus in self: self[locus].calc_consensus()
    
    def getSeqName(self,line):
        firstquote = line.find("'")+1
        name = line[firstquote:line.find("'",firstquote)]
        length = line[line.find("(")+1:line.find(" bp)")]
        return name,length

'''
        print "Working on",self.name
        seqFile = open("tmp/%s_genomicSequence.fasta" % (self.name),"w")
        seqFile.write(">%s\n%s\n" % (self.name,seq))
        seqFile.close()
        newFasta = open("tmp/%s_orfs.fasta" %(self.name),"w")
        orfs = getOrfs("tmp/%s_genomicSequence.fasta" % (self.name),90)
        for seq in orfs.Seqlist: writeFasta(seq, newFasta, "fasta")
        newFasta.close()
        allBlast = open("AnnotateCMDs.sh","a")
        cmd = "bsub -P gis-jbrowse-blastall -J CasBLAST -q prod -R \"rusage[mem=1048,scr=1000]\" " \
        "-R \"span[hosts=1]\" -n 12 -e tmp/%s_blast.err -o tmp/%s_blast.out \"bash tmp/%s.sh; rm -f tmp/%s.sh tmp/%s_genomicSequence.fasta tmp/%s_orfs.fasta\""% (self.name,self.name,self.name,self.name,self.name,self.name)
        blastCMDs = "blastp -db blast/db/nr "\
        "-query tmp/%s_orfs.fasta -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames stitle\" " \
        "-num_threads 12 -out data/Annotations/%s_orfs.blastout\n" % (self.name,self.name)
        cmdFile = open("tmp/%s.sh"%(self.name),"w")
        cmdFile.write(blastCMDs)
        cmdFile.close()
        allBlast.write(cmd+"\n")
        allBlast.close()
'''